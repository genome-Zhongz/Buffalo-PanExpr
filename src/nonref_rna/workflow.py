import concurrent.futures
import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
from .config import resolve_reference,resolve_tools,build_resources
from .io import signature,atomic_json,digest,read_table
from .samples import load_samples,save_samples
from .evidence import summarize,render_report


def check_identity(work,state):
    """Resource-only changes are safe; biological input changes require a new run."""
    work=Path(work);work.mkdir(parents=True,exist_ok=True)
    stable={k:v for k,v in state.items() if k!='resources'}
    stable['inputs']=[signature(p) for p in stable['inputs']]
    fingerprint=hashlib.sha256(json.dumps(stable,sort_keys=True).encode()).hexdigest()
    target=work/'run_manifest.json'
    if target.exists() and json.loads(target.read_text(encoding='utf-8'))['fingerprint']!=fingerprint:
        raise ValueError('Inputs, reference, tools or biological parameters changed; use a NEW output directory. Existing results were not overwritten.')
    atomic_json(target,dict(fingerprint=fingerprint,identity=stable,resources=state.get('resources',{}),timestamp=datetime.datetime.now().astimezone().isoformat()))


def reference_contigs(engine,reference):
    fa=Path(engine['extended_fa']);index=Path(str(fa)+'.fai');lengths={}
    if index.is_file() and index.stat().st_mtime_ns>=fa.stat().st_mtime_ns:
        with index.open() as f:
            for line in f:
                c=line.split('\t');lengths[c[0]]=int(c[1])
    else:
        ident=None
        with fa.open() as f:
            for line in f:
                if line.startswith('>'):
                    ident=line[1:].split()[0]
                    if ident in lengths:raise ValueError('Duplicate reference FASTA ID: '+ident)
                    lengths[ident]=0
                elif line.strip():
                    if ident is None:raise ValueError('Sequence before FASTA header')
                    lengths[ident]+=len(line.strip())
    if reference.get('nonref_contigs'):
        nr=set(Path(reference['nonref_contigs']).read_text().splitlines())
    else:
        prefix=reference.get('nonref_prefix','NR_')
        if not prefix:raise ValueError('Empty nonref_prefix is unsafe; supply an explicit contig list')
        nr={x for x in lengths if x.startswith(prefix)}
    if not nr or nr-set(lengths):raise ValueError('Nonreference contig list/prefix matches zero or unknown contigs')
    for key in ['base_gtf','guide_gtf','nonref_guide']:
        with open(engine[key]) as f:
            for line in f:
                if not line.strip() or line.startswith('#'):continue
                c=line.rstrip().split('\t')
                if len(c)!=9 or c[0] not in lengths:raise ValueError('GTF/FASTA contig disagreement in '+key)
                if int(c[3])<1 or int(c[4])<int(c[3]) or int(c[4])>lengths[c[0]]:raise ValueError('GTF coordinate outside sequence: '+key)
                if key=='base_gtf' and c[0] in nr:raise ValueError('Nonreference contig overlaps original annotation')
                if key=='nonref_guide' and c[0] not in nr:raise ValueError('Prediction guide includes original reference contig')
    return sorted(nr)


def prepare_engine_config(c,work):
    work=Path(work).resolve();r=resolve_reference(c);tools=resolve_tools(c)
    allocation=build_resources(c['resources']['threads'],c['resources']['jobs'])
    return dict(r,tools=tools,work=str(work),**allocation,minimum_free_gib=c['resources']['minimum_free_gib'],
                sort_memory=c['resources']['sort_memory'],compare_original=c['compare_original'],
                pilot_pairs=c['pilot_pairs'],strand_min_reads=c['strand_min_reads'],
                nonref_contigs=str(work/'annotation/nonref_contigs.txt'))


def validate(c,samples_path,work):
    rows=load_samples(samples_path);engine=prepare_engine_config(c,work)
    nr=reference_contigs(engine,c['reference'])
    if c['expression']['min_samples']>len(rows):raise ValueError('expression.min_samples exceeds libraries')
    if c['homology']['enabled'] and not Path(c['homology'].get('database','')).is_file():raise ValueError('homology.enabled requires a local DIAMOND database')
    validate_analysis(rows,c['analysis'])
    return rows,engine,nr


def validate_analysis(rows,options):
    if not options.get('differential'):return
    control=options.get('control');treatment=options.get('treatment')
    if not control or not treatment or control==treatment:raise ValueError('Set distinct analysis.control and analysis.treatment group labels')
    chosen=[r for r in rows if r.get('group') in (control,treatment)]
    for group in [control,treatment]:
        if sum(r.get('group')==group for r in chosen)<2:raise ValueError('Differential analysis requires >=2 independent libraries per group: '+group)
    ids=[r.get('individual') for r in chosen if r.get('individual') not in (None,'','unspecified','unconfirmed','not_provided')]
    if len(ids)!=len(set(ids)):raise ValueError('Repeated individual IDs are incompatible with this unpaired two-group model')
    if not 1<=int(options.get('min_samples',2))<=len(chosen):raise ValueError('analysis.min_samples exceeds selected libraries')


def prepare_reads(c,rows,engine):
    if c['reads']['cleaned']:return rows
    root=Path(engine['work']);dest=root/'qc';dest.mkdir(exist_ok=True)
    def one(row):
        row=dict(row);sid=row['sample'];r1=dest/(sid+'_R1.fastq.gz');r2=dest/(sid+'_R2.fastq.gz');marker=dest/(sid+'.json')
        valid=False
        if marker.exists() and r1.is_file() and r2.is_file():
            saved=json.loads(marker.read_text(encoding='utf-8'))
            valid=saved['outputs']==[signature(r1),signature(r2)]
        if not valid:
            cmd=[engine['tools']['fastp'],'-i',row['r1'],'-I',row['r2'],'-o',str(r1),'-O',str(r2),'-w',str(min(16,engine['stringtie_threads'])),'--json',str(dest/(sid+'.fastp.json')),'--html',str(dest/(sid+'.fastp.html'))]
            with (dest/(sid+'.log')).open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
            atomic_json(marker,dict(command=cmd,outputs=[signature(r1),signature(r2)]))
        row.update(r1=str(r1),r2=str(r2),r1_bytes=r1.stat().st_size,r2_bytes=r2.stat().st_size)
        return row
    with concurrent.futures.ThreadPoolExecutor(max_workers=engine['concurrent_samples']) as pool:return list(pool.map(one,rows))


def preflight_engine(runner):
    from .engine.pipeline import log
    c=runner.c
    for row in runner.samples:
        for key in ['r1','r2']:
            if Path(row[key]).stat().st_size!=int(row[key+'_bytes']):raise ValueError('FASTQ changed: '+row[key])
    for key in ['extended_fa','base_fa']:
        fa=c[key];idx=Path(fa+'.fai')
        if not idx.is_file() or idx.stat().st_mtime_ns<Path(fa).stat().st_mtime_ns:
            runner.run([c['tools']['samtools'],'faidx',fa],runner.root/'logs/faidx.log')
    atomic_json(runner.root/'run_config.snapshot.json',c)
    with (runner.root/'logs/versions.txt').open('w') as f:
        for name,arg in [('hisat2','--version'),('samtools','--version'),('stringtie','--version'),('featureCounts','-v'),('gffcompare','--version'),('gffread','--version')]:
            proc=subprocess.run([c['tools'][name],arg],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,errors='replace')
            f.write(name+': '+proc.stdout+'\n')
    atomic_json(runner.root/'status.json',dict(state='RUNNING',stage='preflight',samples=len(runner.samples)))
    log('Preflight passed: '+str(len(runner.samples))+' paired libraries')


def postprocess(c,work,output=None):
    validate_analysis(read_table(Path(work)/'samples.tsv'),c.get('analysis',{}))
    output=Path(output or work).resolve();report=summarize(work,output,c)
    if c.get('homology',{}).get('enabled'):
        from .homology import annotate
        tools=resolve_tools(c)
        report['homology']=annotate(output/'candidates/longest_transcript_per_gene.fa',c['homology']['database'],tools['diamond'],output/'homology',c['resources']['threads'],c['homology'])
        report['homology_status']='completed'
        candidates=read_table(output/'candidates/nonreference_candidates.tsv')
        byid={r['gene_id']:r for r in read_table(output/'homology/gene_homology.tsv')}
        from .io import write_table
        for r in candidates:r['homology_status']=byid[r['gene_id']]['evidence_class']
        fields=list(candidates[0]) if candidates else list(read_table(output/'candidates/nonreference_candidates.tsv'))
        if candidates:write_table(output/'candidates/nonreference_candidates.tsv',fields,candidates)
    if c.get('analysis',{}).get('plots') or c.get('analysis',{}).get('differential'):
        script=Path(__file__).parent/'scripts/analysis.R'
        params=output/'report/analysis_config.json';atomic_json(params,c['analysis'])
        atomic_json(output/'report/analysis_manifest.json',dict(script_sha256=digest(script),options=c['analysis'],TPM=signature(Path(work)/'results/gene_TPM.tsv')))
        cmd=[resolve_tools(c)['rscript'],str(script),str(Path(work).resolve()),str(output),str(params)]
        with (output/'report/R_analysis.log').open('w') as f:subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT,check=True)
        if c['analysis'].get('differential'):
            from .primary import export_primary_results
            report['primary_differential']=export_primary_results(work,output,c)
    atomic_json(output/'report/summary.json',report);render_report(output,report)
    return report


def run(c,samples_path,work,dry_run=False):
    rows,engine,nr=validate(c,samples_path,work)
    if dry_run:
        return dict(dry_run=True,libraries=len(rows),nonreference_contigs=len(nr),engine_config=engine,
                    steps=['optional_fastp','index','mapping_and_observed_assembly','merge_nonreference','fixed_model_quantification','candidate_evidence_and_FASTA','optional_homology_and_R_analysis'])
    if sys.platform!='linux':raise ValueError('Full workflow requires Linux; validation and summarize work on other platforms')
    work=Path(work).resolve();work.mkdir(parents=True,exist_ok=True)
    import fcntl
    with (work/'pipeline.lock').open('w') as lock:
        try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        except BlockingIOError:raise ValueError('Another run is using this output directory')
        started=False
        try:
            refs=[engine[k] for k in ['extended_fa','guide_gtf','base_gtf','nonref_guide','base_fa']]
            code={str(p.relative_to(Path(__file__).parent)):digest(p) for p in Path(__file__).parent.rglob('*.py')}
            check_identity(work,dict(inputs=sorted(set(refs+[r[k] for r in rows for k in ['r1','r2']]+list(engine['tools'].values()))),
                                    biology=dict(samples=[{k:v for k,v in r.items() if k not in ['r1_bytes','r2_bytes']} for r in rows],cleaned=c['reads']['cleaned'],
                                                 pilot_pairs=c['pilot_pairs'],strand_min_reads=c['strand_min_reads'],compare_original=c['compare_original'],nonref_ids=nr,
                                                 expression=c['expression'],homology=c['homology'],analysis=c['analysis']),
                                    code=code,resources=c['resources']))
            started=True
            (work/'annotation').mkdir(exist_ok=True)
            Path(engine['nonref_contigs']).write_text(''.join(x+'\n' for x in nr))
            prepared=prepare_reads(c,rows,engine);save_samples(work/'samples.tsv',prepared)
            from .engine.pipeline import Pipeline
            runner=Pipeline(engine);runner.run_all()
            result=postprocess(c,work)
            atomic_json(work/'status.json',dict(state='COMPLETED',stage='report',summary=result))
            return result
        except BaseException as e:
            if started:
                (work/'COMPLETED').unlink(missing_ok=True)
                atomic_json(work/'status.json',dict(state='FAILED',error=str(e)))
            raise
