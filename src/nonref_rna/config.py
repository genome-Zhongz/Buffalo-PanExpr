import copy
import json
from pathlib import Path
import shutil
import yaml

DEFAULTS=dict(reference=dict(bundle='',nonref_prefix='NR_'),reads=dict(cleaned=True),
              resources=dict(threads=16,jobs=1,minimum_free_gib=30,sort_memory='768M'),
              expression=dict(tpm=1.0,min_samples=1),homology=dict(enabled=False,evalue=1e-5,max_targets=25),
              analysis=dict(plots=False,differential=False,control='',treatment='',cpm=1,min_samples=2),
              compare_original=False,pilot_pairs=200000,strand_min_reads=1000,tools={})
TOOLS={'hisat2':'hisat2','hisat2_build':'hisat2-build','extract_splice':'hisat2_extract_splice_sites.py',
       'extract_exons':'hisat2_extract_exons.py','samtools':'samtools','stringtie':'stringtie',
       'featureCounts':'featureCounts','gffcompare':'gffcompare','gffread':'gffread'}


def build_resources(threads,jobs):
    if not isinstance(threads,int) or not isinstance(jobs,int) or jobs<1 or threads<4*jobs:
        raise ValueError('Require threads >= 4 * jobs, with positive integer threads/jobs')
    per=threads//jobs;sort=min(2,max(1,per//4));align=per-sort-1
    return dict(index_threads=threads,alignment_threads=align,sort_threads=sort,
                stringtie_threads=per-1,concurrent_samples=jobs)


def load_config(path):
    path=Path(path).resolve()
    with path.open(encoding='utf-8') as f:data=yaml.safe_load(f) or {}
    if not isinstance(data,dict):raise ValueError('Configuration must be a mapping')
    c=copy.deepcopy(DEFAULTS)
    for k,v in data.items():
        if k not in c:raise ValueError('Unknown configuration key: '+k)
        if isinstance(c[k],dict):
            if not isinstance(v,dict):raise ValueError(k+' must be a mapping')
            c[k].update(v)
        else:c[k]=v
    c['_config_path']=str(path)
    for section,keys in [('reference',['bundle','extended_fasta','guide_gtf','original_fasta','original_gtf','nonref_gtf','nonref_contigs']),('homology',['database'])]:
        for key in keys:
            if c[section].get(key):
                p=Path(c[section][key]).expanduser()
                c[section][key]=str((path.parent/p if not p.is_absolute() else p).resolve())
    build_resources(c['resources']['threads'],c['resources']['jobs'])
    e=c['expression']
    if not isinstance(e['min_samples'],int) or e['min_samples']<1 or not 0<=float(e['tpm'])<float('inf'):
        raise ValueError('expression requires finite nonnegative tpm and min_samples >=1')
    if not isinstance(c['reads']['cleaned'],bool):raise ValueError('reads.cleaned must be true or false')
    return c


def resolve_reference(c):
    r=c['reference'];core=Path(r['bundle'])
    if not core.is_dir():raise ValueError('Reference bundle not found: '+str(core))
    old=json.loads((core/'run_config.json').read_text()) if (core/'run_config.json').is_file() else {}
    paths=dict(core=str(core),extended_fa=r.get('extended_fasta',str(core/'rnaseq_reference/extended.fa')),
               guide_gtf=r.get('guide_gtf',str(core/'rnaseq_reference/guide.gtf')),
               base_gtf=r.get('original_gtf',str(core/'reference/reference.normalized.gtf')),
               nonref_guide=r.get('nonref_gtf',str(core/'rnaseq_reference/nonref.guide.gtf')),
               base_fa=r.get('original_fasta',old.get('graph_reference',old.get('annotation_reference',''))))
    for k in ['extended_fa','guide_gtf','base_gtf','nonref_guide']:
        if not Path(paths[k]).is_file():raise ValueError('Missing reference file '+k+': '+paths[k])
    if c['compare_original'] and not Path(paths['base_fa']).is_file():raise ValueError('Missing original_fasta for compare_original')
    if not paths['base_fa'] or not Path(paths['base_fa']).is_file():paths['base_fa']=paths['extended_fa']
    return paths


def resolve_tools(c):
    names=dict(TOOLS)
    if not c['reads']['cleaned']:names['fastp']='fastp'
    if c['homology']['enabled']:names['diamond']='diamond'
    if c['analysis']['plots'] or c['analysis']['differential']:names['rscript']='Rscript'
    result={}
    for key,default in names.items():
        configured=c['tools'].get(key,default);found=shutil.which(configured)
        if not found and Path(configured).is_file():found=str(Path(configured).resolve())
        if not found:raise ValueError('Executable not found: '+key+' ('+configured+'); set tools.'+key+' or PATH')
        result[key]=found
    return result
