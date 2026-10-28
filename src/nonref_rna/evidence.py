"""Observed RNA candidates; explicit expression definitions and sequence export."""
from collections import defaultdict
import html
import math
from pathlib import Path
import re
import shutil
from .io import read_table,write_table,write_fasta,fasta,atomic_json,digest


def union_length(intervals):
    total=0;end=-1
    for a,b in sorted(set(intervals)):
        total+=max(0,b-max(a,end+1)+1);end=max(end,b)
    return total


def summarize(work,output,config,gene_ids=None,copy_matrices=True):
    work=Path(work).resolve();output=Path(output).resolve()
    samples=read_table(work/'samples.tsv');sampleids=[r['sample'] for r in samples]
    if len(set(sampleids))!=len(sampleids) or not samples:raise ValueError('Duplicate/empty samples in result workspace')
    metadata_rows=read_table(work/'results/gene_metadata.tsv')
    metadata={r['gene_id']:r for r in metadata_rows}
    explicit=set(gene_ids) if gene_ids is not None else None
    if explicit is not None:
        for g in explicit:
            if g not in metadata or str(metadata[g]['is_nonreference']).lower() not in ('1','true'):
                raise ValueError('Explicit gene list contains unknown/original-reference gene: '+g)
    if len(metadata)!=len(metadata_rows):raise ValueError('Duplicate metadata gene IDs')
    matrix=read_table(work/'results/gene_TPM.tsv')
    if len({r['gene_id'] for r in matrix})!=len(matrix):raise ValueError('Duplicate TPM gene IDs')
    if set(metadata)!=set(r['gene_id'] for r in matrix):raise ValueError('TPM and annotation genes disagree')
    settings=config.get('expression',{});threshold=float(settings.get('tpm',1));minimum=int(settings.get('min_samples',1))
    if minimum<1 or not math.isfinite(threshold) or threshold<0:raise ValueError('Invalid expression rule')
    if minimum>len(samples):raise ValueError('expression.min_samples exceeds number of libraries')
    tissues=defaultdict(list)
    for row in samples:tissues[row.get('tissue') or 'unspecified'].append(row['sample'])
    membership=defaultdict(set);nrmembership=defaultdict(set)
    candidates=[];tissue_rows=[]
    for row in matrix:
        if set(row)-{'gene_id'}!=set(sampleids):raise ValueError('TPM columns disagree with samples.tsv')
        gid=row['gene_id'];meta=metadata[gid]
        values={s:float(row[s]) for s in sampleids}
        if any(not math.isfinite(v) or v<0 for v in values.values()):raise ValueError('Invalid TPM for '+gid)
        nr=str(meta['is_nonreference']).lower() in ('1','true')
        support=sum(v>=threshold for v in values.values())
        for tissue,ids in tissues.items():
            mean=sum(values[s] for s in ids)/len(ids)
            if mean>1:
                membership[tissue].add(gid)
                if nr:nrmembership[tissue].add(gid)
                tissue_rows.append(dict(gene_id=gid,tissue=tissue,mean_TPM=mean,is_nonreference=int(nr)))
        if nr and (support>=minimum if explicit is None else gid in explicit):
            candidates.append(dict(meta,samples_TPM_ge_threshold=support,max_TPM=max(values.values()),mean_TPM=sum(values.values())/len(values)))
    wanted={r['gene_id'] for r in candidates};models={};exons=defaultdict(list);gtflines=[];nr_contigs=set()
    with (work/'annotation/nonref.rna_merged.gtf').open() as f:
        for line in f:
            if line.startswith('#') or not line.strip():continue
            c=line.rstrip().split('\t')
            if len(c)!=9:raise ValueError('Malformed nonreference GTF')
            nr_contigs.add(c[0])
            a=dict(re.findall(r'(\w+) "([^"]*)"',c[8]));g=a.get('gene_id');t=a.get('transcript_id')
            if g not in wanted:continue
            if not t:raise ValueError('Missing transcript ID: '+g)
            gtflines.append(line)
            if c[2]=='transcript':
                if t in models:raise ValueError('Duplicate transcript '+t)
                models[t]=dict(gene_id=g,transcript_id=t,contig=c[0],strand=c[6],start=int(c[3]),end=int(c[4]))
            if c[2]=='exon':exons[t].append((int(c[3]),int(c[4])))
    seqs={}
    for title,seq in fasta(work/'annotation/nonref.rna_transcripts.fa'):
        tid=title.split()[0]
        if tid not in models:continue
        if tid in seqs:raise ValueError('Duplicate FASTA transcript '+tid)
        seqs[tid]=seq
    if set(seqs)!=set(models):raise ValueError('Missing transcript sequence(s): '+','.join(sorted(set(models)-set(seqs))[:10]))
    bygene=defaultdict(list)
    for tid,m in models.items():
        if not seqs[tid] or len(seqs[tid])!=union_length(exons[tid]):raise ValueError('Sequence/exon length disagreement: '+tid)
        if set(seqs[tid].upper())-set('ACGTRYSWKMBDHVN'):raise ValueError('Invalid nucleotide alphabet: '+tid)
        m.update(length_bp=len(seqs[tid]),N_percent=100*seqs[tid].upper().count('N')/len(seqs[tid]))
        bygene[m['gene_id']].append(tid)
    if set(bygene)!=wanted:raise ValueError('Candidate gene missing observed transcript: '+','.join(wanted-set(bygene)))
    from .provenance import verify_annotation_and_quant,verify_transcript_sequences
    verify_annotation_and_quant(work,metadata,matrix,samples,nr_contigs,models,exons)
    source_reference=verify_transcript_sequences(work,models,exons,seqs)
    longest={g:sorted(ts,key=lambda t:(-len(seqs[t]),t))[0] for g,ts in bygene.items()}
    structure=defaultdict(list)
    for path in sorted((work/'annotation').glob('nr_vs_augustus*.tmap')):
        for r in read_table(path):structure[r['qry_gene_id']].append(r)
    for row in candidates:
        g=row['gene_id'];tid=longest[g];hits=structure[g]
        row.update(representative_transcript=tid,transcript_count=len(bygene[g]),representative_length_bp=len(seqs[tid]),
                   N_percent=models[tid]['N_percent'],gffcompare_codes=';'.join(sorted({h['class_code'] for h in hits})),
                   augustus_gene_ids=';'.join(sorted({h['ref_gene_id'] for h in hits if h['ref_gene_id'] not in ('-','.','')})),
                   augustus_structure_equal=int(any(h['class_code']=='=' for h in hits)),
                   homology_status='not_run',interpretation='RNA-supported nonreference candidate; coding function not established')
    dest=output/'candidates';dest.mkdir(parents=True,exist_ok=True)
    fields=['gene_id','contig','strand','start','end','gene_name','samples_TPM_ge_threshold','max_TPM','mean_TPM','transcript_count','representative_transcript','representative_length_bp','N_percent','gffcompare_codes','augustus_gene_ids','augustus_structure_equal','homology_status','interpretation']
    write_table(dest/'nonreference_candidates.tsv',fields,candidates)
    with (dest/'all_transcripts.fa').open('w') as allfa,(dest/'longest_transcript_per_gene.fa').open('w') as repfa:
        for row in candidates:
            g=row['gene_id']
            for tid in sorted(bygene[g]):write_fasta(allfa,f'{tid} gene_id={g} length_bp={len(seqs[tid])}',seqs[tid])
            tid=longest[g];write_fasta(repfa,f'{g} transcript_id={tid} length_bp={len(seqs[tid])}',seqs[tid])
    (dest/'selected_models.gtf').write_text(''.join(gtflines))
    (dest/'gene_ids.txt').write_text(''.join(r['gene_id']+'\n' for r in candidates))
    for tid,m in models.items():m['representative_longest']=int(longest[m['gene_id']]==tid)
    write_table(dest/'gene_transcript_mapping.tsv',['gene_id','transcript_id','contig','strand','start','end','length_bp','N_percent','representative_longest'],models.values())
    write_table(output/'report/genes_tissue_mean_TPM_gt1.tsv',['gene_id','tissue','mean_TPM','is_nonreference'],tissue_rows)
    write_table(output/'report/tissue_expression_counts.tsv',['tissue','all_genes','nonreference_genes'],[dict(tissue=t,all_genes=len(membership[t]),nonreference_genes=len(nrmembership[t])) for t in tissues])
    if copy_matrices:
        matrices=output/'matrices';matrices.mkdir(exist_ok=True)
        for name in ['gene_TPM.tsv','nonref_TPM.tsv','gene_counts.raw.tsv','nonref_counts.raw.tsv']:
            src=work/'results'/name
            if src.is_file() and src.resolve()!=(matrices/name).resolve():shutil.copyfile(src,matrices/name)
    report=dict(libraries=len(samples),annotated_genes=len(matrix),candidate_genes=len(candidates),candidate_transcripts=len(models),
                expression_rule=f'TPM >= {threshold}; at least {minimum} libraries; unique genes' if explicit is None else 'Explicit differential candidate list; TPM support reported but not used to filter again',
                tissue_rule='Arithmetic tissue mean TPM > 1; union deduplicated across tissues',
                all_genes_tissue_mean_gt1_union=len(set().union(*membership.values())),
                nonref_genes_tissue_mean_gt1_union=len(set().union(*nrmembership.values())),
                nonref_tissue_sum_with_duplicates=sum(map(len,nrmembership.values())),
                candidates_with_augustus_structure_equal=sum(r['augustus_structure_equal'] for r in candidates),
                homology_status='not_run',source_work=str(work),
                source_matrix_sha256=digest(work/'results/gene_TPM.tsv'),
                source_GTF_sha256=digest(work/'annotation/nonref.rna_merged.gtf'),
                source_transcript_FASTA_sha256=digest(work/'annotation/nonref.rna_transcripts.fa'),
                source_reference=source_reference,source_provenance='Counting GTF, source per-sample quantification, reference manifest and transcript/reference sequences verified',
                limitations=['Expression and assembly do not prove protein-coding function.',
                             'gffcompare classes are relative to the supplied prediction guide, not novelty against the original genome.',
                             'Homologs may represent related copies, fragments or pseudogenes; these states require independent validation.'])
    atomic_json(output/'report/summary.json',report)
    render_report(output,report)
    return report


def render_report(output,summary):
    output=Path(output);rows=read_table(output/'candidates/nonreference_candidates.tsv')
    keys=['gene_id','contig','samples_TPM_ge_threshold','max_TPM','representative_transcript','augustus_structure_equal']
    cells='<table><tr>'+''.join('<th>'+k+'</th>' for k in keys)+'</tr>'
    for r in rows:cells+='<tr>'+''.join('<td>'+html.escape(str(r.get(k,'')))+'</td>' for k in keys)+'</tr>'
    content='<!doctype html><meta charset="utf-8"><title>nonref-rna results</title><style>body{font:15px sans-serif;margin:2em}table{border-collapse:collapse}td,th{border:1px solid #ddd;padding:6px}th{background:#eee}</style><h1>Nonreference RNA candidates</h1><pre>'+html.escape(__import__('json').dumps(summary,indent=2,ensure_ascii=False))+'</pre>'+cells+'</table>'
    if summary.get('primary_differential'):
        content=content.replace('<h1>Nonreference RNA candidates</h1>','<h1>Nonreference RNA candidates</h1><p><strong>Primary differential results: raw P &lt;0.05 and |log2FC| &gt;1 (exploratory).</strong> <a href="../primary_Pvalue/nonreference_candidates_P05.csv">Open P-value candidate CSV</a> · <a href="../primary_Pvalue/nonreference_candidates_P05.longest.fa">Download candidate FASTA</a></p>')
    for label,stats in [('Expression candidates',summary.get('homology')),('Raw-P differential candidates',summary.get('primary_differential',{}).get('homology'))]:
        if stats and 'text_zh' in stats:
            content+='<h2>'+label+'</h2><p>'+html.escape(stats['text_zh'])+'</p><p>'+html.escape(stats['text_en'])+'</p>'
    (output/'report/report.html').write_text(content,encoding='utf-8')
