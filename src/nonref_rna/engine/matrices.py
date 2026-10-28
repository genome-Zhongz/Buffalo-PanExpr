#!/usr/bin/env python3
import collections
import csv
import json
from pathlib import Path
import re


def parse_abundance(path):
    values=collections.defaultdict(float)
    with open(path) as inp:
        reader=csv.DictReader(inp,delimiter='\t')
        if not {'Gene ID','TPM'}<=set(reader.fieldnames or []):raise ValueError('Unexpected StringTie abundance header')
        for row in reader:
            gid=row['Gene ID']
            if gid and gid not in ('-','.'):values[gid]+=float(row['TPM'])
    return dict(values)


def parse_featurecounts(path):
    # Highly alternatively spliced reference genes can exceed csv's 128 KiB default.
    csv.field_size_limit(100000000)
    counts,info={},{}
    with open(path) as inp:
        reader=csv.DictReader((line for line in inp if not line.startswith('#')),delimiter='\t')
        fields=reader.fieldnames or []
        if len(fields)!=7 or fields[0]!='Geneid':raise ValueError('Expected single-sample featureCounts table')
        for row in reader:
            gid=row['Geneid']
            if gid in counts:raise ValueError('Duplicate gene in count table: '+gid)
            value=float(row[fields[-1]])
            if value<0 or value!=int(value):raise ValueError('Noninteger count')
            counts[gid]=int(value)
            info[gid]={k:row[k] for k in fields[1:6]}
    return counts,info


def write_matrix(path,samples,genes,values):
    with open(path,'w',newline='') as out:
        writer=csv.writer(out,delimiter='\t',lineterminator='\n')
        writer.writerow(['gene_id']+samples)
        for gid in genes:
            row=[gid]
            for sample in samples:
                value=values[sample].get(gid,0)
                row.append(format(value,'.8g') if isinstance(value,float) else str(value))
            writer.writerow(row)


def select_nonref_gtf(source,nrids,destination):
    with open(source) as inp,open(destination,'w') as out:
        for line in inp:
            if line.startswith('#') or not line.strip():continue
            if line.split('\t',1)[0] in nrids:out.write(line)


def gtf_genes(path,nonref_ids=None):
    genes={}
    with open(path) as inp:
        for line in inp:
            if line.startswith('#') or not line.strip():continue
            c=line.rstrip('\n').split('\t')
            attrs=dict(re.findall(r'(\w+) "([^"]*)"',c[8]))
            gid=attrs.get('gene_id')
            if not gid:raise ValueError('Missing gene_id in final annotation')
            location=(c[0],c[6])
            if gid in genes and (genes[gid]['contig'],genes[gid]['strand'])!=location:raise ValueError('Gene ID collision: '+gid)
            m=genes.setdefault(gid,dict(contig=c[0],strand=c[6],start=int(c[3]),end=int(c[4]),gene_name=attrs.get('gene_name',''),is_nonreference=(c[0] in nonref_ids if nonref_ids is not None else c[0].startswith('NR_'))))
            m['start']=min(m['start'],int(c[3]));m['end']=max(m['end'],int(c[4]))
    return genes


def collect(work,samples,core,nonref_ids=None):
    work=Path(work);dest=work/'results';dest.mkdir(exist_ok=True)
    genes=gtf_genes(work/'annotation/final_annotation.gtf',nonref_ids);geneids=sorted(genes);sampleids=[r['sample'] for r in samples]
    tpm,counts,info={}, {},None
    for row in samples:
        sid=row['sample']
        tpm[sid]=parse_abundance(work/'quant'/sid/'gene_abundance.tsv')
        counts[sid],this_info=parse_featurecounts(work/'counts'/f'{sid}.tsv')
        if set(counts[sid])!=set(genes):raise ValueError('Count genes disagree with annotation: '+sid)
        if info is not None and this_info!=info:raise ValueError('Count gene coordinates differ: '+sid)
        info=this_info
        if sum(counts[sid].values())==0:raise ValueError('Zero assigned fragments: '+sid)
        extra=set(tpm[sid])-set(genes)
        if extra:raise ValueError('Unexpected TPM gene IDs: '+sid)
    nrgenes=[g for g in geneids if genes[g]['is_nonreference']]
    write_matrix(dest/'gene_TPM.tsv',sampleids,geneids,tpm)
    write_matrix(dest/'nonref_TPM.tsv',sampleids,nrgenes,tpm)
    write_matrix(dest/'gene_counts.raw.tsv',sampleids,geneids,counts)
    write_matrix(dest/'nonref_counts.raw.tsv',sampleids,nrgenes,counts)
    totals={s:sum(counts[s].values()) for s in sampleids}
    cpm={s:{g:v*1e6/totals[s] for g,v in counts[s].items()} for s in sampleids}
    write_matrix(dest/'nonref_CPM.descriptive.tsv',sampleids,nrgenes,cpm)
    with open(dest/'gene_metadata.tsv','w',newline='') as out:
        writer=csv.writer(out,delimiter='\t');writer.writerow(['gene_id','contig','strand','start','end','gene_name','is_nonreference'])
        for gid in geneids:
            m=genes[gid];writer.writerow([gid,m['contig'],m['strand'],m['start'],m['end'],m['gene_name'],int(m['is_nonreference'])])
    mapping=collections.defaultdict(list)
    for path in (work/'annotation').glob('nr_vs_augustus*.tmap'):
        with open(path) as inp:
            for row in csv.DictReader(inp,delimiter='\t'):
                mapping[row['qry_gene_id']].append((row['ref_gene_id'],row['ref_id'],row['class_code']))
    support_file=Path(core)/'blastp/gene_ids.80_80.txt'
    blastp=set(support_file.read_text().splitlines()) if support_file.is_file() else set()
    groups=collections.defaultdict(list)
    for row in samples:
        groups['tissue:'+row['tissue']].append(row['sample'])
        groups['tissue_breed:'+row['tissue']+'/'+row['breed']].append(row['sample'])
    means={group:{g:sum(tpm[s].get(g,0) for s in members)/len(members) for g in nrgenes} for group,members in groups.items()}
    write_matrix(dest/'nonref_TPM.group_means.tsv',list(groups),nrgenes,means)
    with open(dest/'nonref_evidence.tsv','w',newline='') as out:
        writer=csv.writer(out,delimiter='\t')
        writer.writerow(['gene_id','contig','samples_TPM_ge_1','max_sample_TPM','tissues_mean_TPM_ge_1','tissue_breed_mean_TPM_ge_1','augustus_gene_ids','gffcompare_codes','structural_match_equal','related_structure_c_or_j','blastp_supported_prediction','expressed_and_structural_match'])
        strict=0
        for gid in nrgenes:
            hits=mapping[gid];refs={r for r,tid,c in hits if r not in ('-','.','')};codes={c for r,tid,c in hits}
            expressed=sum(tpm[s].get(gid,0)>=1 for s in sampleids)
            tissues=[group.split(':',1)[1] for group in groups if group.startswith('tissue:') and means[group][gid]>=1]
            breedgroups=[group.split(':',1)[1] for group in groups if group.startswith('tissue_breed:') and means[group][gid]>=1]
            strict+=int(bool(breedgroups) and '=' in codes)
            writer.writerow([gid,genes[gid]['contig'],expressed,max(tpm[s].get(gid,0) for s in sampleids),';'.join(tissues),';'.join(breedgroups),';'.join(sorted(refs)),';'.join(sorted(codes)),int('=' in codes),int(bool(codes&{'c','j'})),int(bool(refs&blastp)),int(bool(breedgroups) and '=' in codes)])
    with open(dest/'sample_QC.tsv','w',newline='') as out:
        writer=csv.writer(out,delimiter='\t');writer.writerow(['sample','tissue','breed','strand','strand_status','forward_fraction','reverse_fraction','informative_read1','extended_alignment_rate_percent','original_alignment_rate_percent','assigned_fragments','NR_genes_TPM_ge_1'])
        for row in samples:
            sid=row['sample'];strand=json.loads((work/'strand'/f'{sid}.json').read_text())
            rates=[]
            for mode in ('extended','original'):
                path=work/'logs'/f'{sid}.{mode}.hisat2.log'
                text=path.read_text() if path.is_file() else ''
                found=re.findall(r'([0-9.]+)% overall alignment rate',text)
                rates.append(found[-1] if found else '')
            writer.writerow([sid,row['tissue'],row['breed'],strand['strand'],strand['status'],strand['forward_fraction'],strand['reverse_fraction'],strand['informative_read1'],*rates,sum(counts[sid].values()),sum(tpm[sid].get(g,0)>=1 for g in nrgenes)])
    summary=dict(samples=len(samples),tissues=dict(collections.Counter(r['tissue'] for r in samples)),genes=len(genes),observed_nonreference_genes=len(nrgenes),
                 nonreference_genes_TPM_ge_1_any_sample=sum(any(tpm[s].get(g,0)>=1 for s in sampleids) for g in nrgenes),
                 nonreference_genes_TPM_ge_1_tissue_breed_mean_and_structural_equal=strict,
                 expression_threshold='TPM >= 1; sample, tissue mean and tissue/breed mean results separately reported',
                 differential_expression_run=False,individual_pairing_confirmed=False)
    (dest/'summary.json').write_text(json.dumps(summary,indent=2))
    return summary
