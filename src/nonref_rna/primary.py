"""P-value primary results; full FDR statistics remain in the DE folder."""
import csv
from pathlib import Path
import shutil
from .io import read_table,atomic_json
from .evidence import summarize
from .homology import summarize_homology,write_homology_summary


def export_primary_results(work,output,config):
    output=Path(output);de=output/'differential_expression';dest=output/'primary_Pvalue'
    selected=read_table(de/'P05_FC1.nonreference_selected.tsv')
    if any(float(r['PValue'])>=.05 or abs(float(r['logFC']))<=1 for r in selected):raise ValueError('Primary P-value table contains an unselected gene')
    selected.sort(key=lambda r:float(r['PValue']))
    genes=[r['gene_id'] for r in selected]
    if len(genes)!=len(set(genes)):raise ValueError('Duplicate primary candidate gene ID')
    summarize(work,dest,config,gene_ids=genes,copy_matrices=False)
    mapping={r['gene_id']:r for r in read_table(dest/'candidates/nonreference_candidates.tsv')}
    homology={}
    if (output/'homology/gene_homology.tsv').is_file():homology={r['gene_id']:r for r in read_table(output/'homology/gene_homology.tsv')}
    fields=['gene_id','log2FC','P','direction','protein_description','matched_accession','homology_evidence','representative_transcript','length_bp']
    with (dest/'nonreference_candidates_P05.csv').open('w',encoding='utf-8-sig',newline='') as f:
        writer=csv.DictWriter(f,fields);writer.writeheader()
        for row in selected:
            g=row['gene_id'];h=homology.get(g,{})
            writer.writerow(dict(gene_id=g,log2FC=row['logFC'],P=row['PValue'],direction=row.get('direction','Up' if float(row['logFC'])>0 else 'Down'),
                                 protein_description=h.get('protein_description') or ('No match in configured database' if h.get('evidence_class')=='no_detected_protein_match' else 'No protein annotation available'),
                                 matched_accession=h.get('matched_accession',''),homology_evidence=h.get('evidence_class','not_run'),representative_transcript=mapping[g]['representative_transcript'],length_bp=mapping[g]['representative_length_bp']))
    full=read_table(de/'P05_FC1.all_selected.tsv')
    with (dest/'all_genes_P05.csv').open('w',encoding='utf-8-sig',newline='') as f:
        fields_all=['gene_id','logFC','P','direction','is_nonreference'];writer=csv.DictWriter(f,fields_all,extrasaction='ignore');writer.writeheader();writer.writerows(dict(r,P=r['PValue']) for r in full)
    aliases={'longest_transcript_per_gene.fa':'nonreference_candidates_P05.longest.fa','all_transcripts.fa':'nonreference_candidates_P05.all_transcripts.fa','gene_transcript_mapping.tsv':'gene_transcript_mapping.tsv','selected_models.gtf':'nonreference_candidates_P05.gtf'}
    for source,name in aliases.items():shutil.copyfile(dest/'candidates'/source,dest/name)
    for source,name in [('Volcano_P05_FC1','Volcano_Pvalue'),('Volcano_P05_FC1_unmarked','Volcano_P_unmarked'),('MDS_nonreference_P_candidates_TPM','MDS_TPM')]:
        for extension in ['png','pdf']:
            p=output/'figures'/(source+'.'+extension)
            if p.is_file():shutil.copyfile(p,dest/(name+'.'+extension))
    report=dict(rule='P < 0.05 and abs(log2FC) > 1',all_gene_candidates=len(full),
                nonreference_candidates=len(genes),nonreference_up=sum(float(r['logFC'])>0 for r in selected),nonreference_down=sum(float(r['logFC'])<0 for r in selected))
    report['homology']=summarize_homology(list(homology.values()),genes)
    write_homology_summary(dest,report['homology'])
    atomic_json(dest/'summary.json',report)
    (dest/'README_zh.txt').write_text('主结果使用P<0.05且|log2FC|>1筛选，按P值排序。\nnonreference_candidates_P05.csv：非参考候选，含P值、方向和同源描述；FASTA提供每基因最长转录本及所有转录本。\nall_genes_P05.csv：全部基因的P值候选。\nVolcano_P_unmarked.png：不带非参考基因标记的火山图。\nFDR统计另存于上级differential_expression目录。\n',encoding='utf-8')
    return report
