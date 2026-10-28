from collections import Counter,defaultdict
import csv
from pathlib import Path
import re
import subprocess
from .io import fasta,write_table,atomic_json

FIELDS='qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore qcovhsp scovhsp qframe stitle'.split()


def summarize_homology(rows, gene_ids):
    """Count one annotation per selected gene, keeping missing annotation explicit."""
    genes=list(gene_ids)
    if len(genes)!=len(set(genes)):
        raise ValueError('Duplicate candidate gene ID in homology summary')
    byid={}
    for row in rows:
        if row['gene_id'] in byid:
            raise ValueError('Duplicate gene ID in homology annotation')
        byid[row['gene_id']]=row
    counts=Counter(byid.get(g,{}).get('evidence_class','not_run') for g in genes)
    broad=counts['high_identity_broad_protein_coverage']
    partial=counts['high_identity_partial_protein_match']
    limited=counts['limited_or_divergent_protein_match']
    known={'high_identity_broad_protein_coverage','high_identity_partial_protein_match',
           'limited_or_divergent_protein_match','no_detected_protein_match','not_run'}
    if set(counts)-known:
        raise ValueError('Unknown homology evidence class')
    n=len(genes);matched=broad+partial+limited
    zh=f'{n} 个候选中，{matched} 个检出蛋白同源匹配，其中 {broad} 个达到较高一致性和较广蛋白覆盖标准。'
    zh+=f'另有 {partial} 个高一致性部分匹配、{limited} 个较短或较低一致性匹配；{counts["no_detected_protein_match"]} 个搜索后未检出匹配，{counts["not_run"]} 个未注释。'
    en=f'Protein homology searches detected matches for {matched} of {n} candidates, including {broad} with high sequence identity and broad protein coverage.'
    en+=f' There were {partial} high-identity partial matches, {limited} limited or divergent matches, {counts["no_detected_protein_match"]} searched candidates without detected matches and {counts["not_run"]} candidates without annotation.'
    return dict(candidates=n,with_matches=matched,high_identity_broad_protein_coverage=broad,
                high_identity_partial_protein_match=partial,limited_or_divergent_protein_match=limited,
                without_matches=counts['no_detected_protein_match'],not_annotated=counts['not_run'],
                text_zh=zh,text_en=en,
                broad_match_rule='identity >=80%; target protein coverage >=80%; alignment >=50 aa; E-value <=1e-10')


def write_homology_summary(output, summary):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    atomic_json(output/'homology_summary.json',summary)
    (output/'homology_summary.md').write_text(
        '# Protein homology summary\n\n'+summary['text_zh']+'\n\n'+summary['text_en']+
        '\n\nBroad-match criteria: '+summary['broad_match_rule']+
        '. Coverage refers to the matched target protein, not transcript completeness.\n'
        '\nCounts describe protein sequence similarity, not AUGUSTUS/assembly dual support or experimentally validated function.\n',encoding='utf-8')


def classify_hit(h):
    if float(h['pident'])>=80 and int(h['length'])>=50 and float(h['evalue'])<=1e-10:
        return 'high_identity_broad_protein_coverage' if float(h['scovhsp'])>=80 else 'high_identity_partial_protein_match'
    return 'limited_or_divergent_protein_match'


def annotate(query,database,diamond,output,threads,options):
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    queries=[title.split()[0] for title,seq in fasta(query)]
    cmd=[diamond,'blastx','--query',str(query),'--db',str(database),'--threads',str(threads),'--ultra-sensitive','--evalue',str(options.get('evalue',1e-5)),
         '--max-target-seqs',str(options.get('max_targets',25)),'--max-hsps','1','--outfmt','6',*FIELDS,'--out',str(output/'blastx.tsv')]
    version=subprocess.run([diamond,'version'],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=True).stdout
    atomic_json(output/'command.json',dict(command=cmd,diamond_version=version,columns=FIELDS))
    if queries:
        with (output/'diamond.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
    else:(output/'blastx.tsv').write_text('')
    hits=defaultdict(list)
    with (output/'blastx.tsv').open() as f:
        for r in csv.reader(f,delimiter='\t'):
            if len(r)!=len(FIELDS):raise ValueError('Malformed BLASTX output')
            h=dict(zip(FIELDS,r));hits[h['qseqid']].append(h)
    rows=[]
    for g in queries:
        hs=sorted(hits[g],key=lambda h:(-float(h['bitscore']),float(h['evalue']),h['sseqid']))
        if not hs:
            rows.append(dict(gene_id=g,evidence_class='no_detected_protein_match',interpretation='No match in configured database; not proof of noncoding status'));continue
        # Keep the score winner explicit; prefer non-PDB accessions only at >=95% of its score.
        pool=[h for h in hs if re.match(r'[NX]P_\d+\.\d+$',h['sseqid']) and float(h['bitscore'])>=.95*float(hs[0]['bitscore'])]
        best=(pool or hs)[0]
        gene=re.search(r'\bGN=([^ ;]+)',best['stitle'])
        rows.append(dict(gene_id=g,matched_accession=best['sseqid'],protein_description=best['stitle'],
                         candidate_symbol=gene[1] if gene else '',identity_percent=best['pident'],query_coverage_percent=best['qcovhsp'],
                         subject_coverage_percent=best['scovhsp'],alignment_aa=best['length'],E_value=best['evalue'],
                         top_scoring_accession=hs[0]['sseqid'],evidence_class=classify_hit(best),
                         interpretation='Homologous protein match; locus-specific orthology and coding function unconfirmed'))
    aliases=Counter(r.get('matched_accession') for r in rows if r.get('matched_accession'))
    for r in rows:r['candidates_sharing_matched_accession']=aliases.get(r.get('matched_accession'),0)
    fields=['gene_id','candidate_symbol','matched_accession','protein_description','identity_percent','query_coverage_percent','subject_coverage_percent','alignment_aa','E_value','top_scoring_accession','candidates_sharing_matched_accession','evidence_class','interpretation']
    write_table(output/'gene_homology.tsv',fields,rows)
    result=dict(queries=len(queries),with_matches=sum(bool(hits[g]) for g in queries),classes=dict(Counter(r['evidence_class'] for r in rows)))
    result.update(summarize_homology(rows,queries))
    write_homology_summary(output,result)
    atomic_json(output/'summary.json',result)
    return result
