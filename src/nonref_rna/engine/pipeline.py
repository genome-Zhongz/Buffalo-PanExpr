#!/usr/bin/env python3
import argparse
import concurrent.futures
from contextlib import contextmanager
import csv
import datetime
import gzip
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import threading
from ..samples import first_record, read_id
from ..io import open_fastq
from .strand import ExonIndex,count_orientation,decide_strand
from .matrices import collect,select_nonref_gtf,gtf_genes


PRINT_LOCK=threading.Lock()


def splice_index_options(splices,exons):
    # HISAT2 2.2.1 raises Nongraph for an empty --ss file; a linear index is
    # appropriate when the supplied annotation contains only single-exon models.
    return ['--ss',str(splices),'--exon',str(exons)] if Path(splices).stat().st_size else []


def log(message):
    with PRINT_LOCK:
        print('['+datetime.datetime.now().astimezone().isoformat(timespec='seconds')+'] '+message,flush=True)


def atomic_json(path,data):
    path=Path(path);temp=path.with_name(path.name+'.tmp.'+str(threading.get_ident()))
    temp.write_text(json.dumps(data,indent=2));os.replace(temp,path)


def file_signature(path):
    path=Path(path);st=path.stat()
    return dict(path=str(path.resolve()),size=st.st_size,mtime_ns=st.st_mtime_ns)


def flags(strand,program):
    if program=='hisat2':return [] if strand=='unstranded' else ['--rna-strandness',{'rf':'RF','fr':'FR'}[strand]]
    if program=='stringtie':return [] if strand=='unstranded' else ['--'+strand]
    if program=='featureCounts':return ['-s',{'unstranded':'0','fr':'1','rf':'2'}[strand]]
    raise ValueError(program)


class Pipeline:
    def __init__(self,config):
        # HISAT2-build uses /usr/bin/env python even when launched by absolute path.
        os.environ['PATH']=str(Path(sys.executable).parent)+os.pathsep+os.environ.get('PATH','')
        self.c=config;self.root=Path(config['work']).resolve();self.t=config['tools'];self.stop=threading.Event()
        self.root.mkdir(exist_ok=True)
        for d in ('index','logs','strand','assembly','alignments','nonref_alignments','scratch','annotation','quant','counts','checkpoints','results','smoke_test'):
            (self.root/d).mkdir(exist_ok=True)
        with open(self.root/'samples.tsv') as inp:self.samples=list(csv.DictReader(inp,delimiter='\t'))
        if not self.samples:raise ValueError('No samples')
        self.nonref_ids=set(Path(config['nonref_contigs']).read_text().splitlines()) if config.get('nonref_contigs') else set()
        self.exons=None;self.disk_condition=threading.Condition();self.disk_reserved=0

    @contextmanager
    def reserve_disk(self,row):
        estimate=max(30*2**30,3*(int(row['r1_bytes'])+int(row['r2_bytes'])))
        with self.disk_condition:
            while True:
                free=shutil.disk_usage(self.root).free;floor=self.c['minimum_free_gib']*2**30
                if free<floor+estimate:raise RuntimeError('Insufficient disk headroom for '+row['sample'])
                if free>=floor+estimate+self.disk_reserved:
                    self.disk_reserved+=estimate;break
                self.disk_condition.wait(timeout=10)
                if self.stop.is_set():raise RuntimeError('Cancelled while waiting for disk reservation')
        try:yield
        finally:
            with self.disk_condition:
                self.disk_reserved-=estimate;self.disk_condition.notify_all()

    def checkpoint_valid(self,marker):
        if not marker.is_file():return False
        record=json.loads(marker.read_text())
        artifacts=record.get('artifacts')
        if not artifacts:return False
        return all(Path(s['path']).is_file() and file_signature(s['path'])==s for s in artifacts)

    def checkpoint(self,marker,paths,**extra):
        atomic_json(marker,dict(artifacts=[file_signature(p) for p in paths],**extra))

    def run(self,cmd,logpath,stdout=None):
        with open(logpath,'a') as err:
            err.write('\nCOMMAND: '+repr([str(x) for x in cmd])+'\n');err.flush()
            if stdout is None:
                subprocess.run([str(x) for x in cmd],stdout=err,stderr=err,check=True)
            else:
                with open(stdout,'w') as out:subprocess.run([str(x) for x in cmd],stdout=out,stderr=err,check=True)

    def pipe(self,first,second,first_log,second_log,stdout=None):
        with open(first_log,'w') as a,open(second_log,'w') as b:
            dest=open(stdout,'w') if stdout else b
            p1=subprocess.Popen([str(x) for x in first],stdout=subprocess.PIPE,stderr=a)
            try:
                p2=subprocess.Popen([str(x) for x in second],stdin=p1.stdout,stdout=dest,stderr=b)
                p1.stdout.close();rc2=p2.wait();rc1=p1.wait()
                if rc1 or rc2:raise RuntimeError(f'Pipeline exit codes {rc1}/{rc2}; logs: {first_log}, {second_log}')
            except BaseException:
                if p1.poll() is None:p1.terminate();p1.wait()
                raise
            finally:
                if stdout:dest.close()

    def safe_unlink(self,path):
        path=Path(path)
        if path.exists():
            if self.root not in path.resolve().parents:raise ValueError('Refusing to delete outside workspace: '+str(path))
            path.unlink(missing_ok=True)

    def space(self,pair_bytes=0):
        needed=max(self.c['minimum_free_gib']*2**30,int(pair_bytes*2.5))
        free=shutil.disk_usage(self.root).free
        if free<needed:raise RuntimeError(f'Insufficient disk: {free/2**30:.1f} GiB available; need {needed/2**30:.1f} GiB')

    def preflight(self):
        from ..workflow import preflight_engine
        preflight_engine(self)

    def build_indexes(self):
        import fcntl
        with open(self.root/'index/build.lock','w') as lock:
            fcntl.flock(lock,fcntl.LOCK_EX)
            for name,fa,gtf in [('extended',self.c['extended_fa'],self.c['guide_gtf']),('original',self.c['base_fa'],self.c['base_gtf'])]:
                if name=='original' and not self.c['compare_original']:continue
                marker=self.root/'checkpoints'/f'index_{name}.json';prefix=self.root/'index'/name
                if self.checkpoint_valid(marker):
                    suffix='ht2' if Path(str(prefix)+'.1.ht2').is_file() else 'ht2l'
                    if not all(Path(str(prefix)+f'.{n}.{suffix}').is_file() for n in range(1,9)):raise ValueError('Incomplete saved index: '+name)
                    continue
                log('Building '+name+' HISAT2 splice-aware index; detailed log: logs/index_'+name+'.log')
                atomic_json(self.root/'status.json',dict(state='RUNNING',stage='index_'+name))
                ss=self.root/'index'/f'{name}.ss';ex=self.root/'index'/f'{name}.exon'
                self.run([sys.executable,self.t['extract_splice'],gtf],self.root/'logs'/f'index_{name}.extract.log',stdout=ss)
                self.run([sys.executable,self.t['extract_exons'],gtf],self.root/'logs'/f'index_{name}.extract.log',stdout=ex)
                self.space()
                self.run([self.t['hisat2_build'],'-p',str(self.c['index_threads']),*splice_index_options(ss,ex),fa,prefix],self.root/'logs'/f'index_{name}.log')
                suffix='ht2' if Path(str(prefix)+'.1.ht2').is_file() else 'ht2l'
                if not all(Path(str(prefix)+f'.{n}.{suffix}').is_file() and Path(str(prefix)+f'.{n}.{suffix}').stat().st_size for n in range(1,9)):
                    raise ValueError('Index build did not generate all files: '+name)
                self.checkpoint(marker,[Path(str(prefix)+f'.{n}.{suffix}') for n in range(1,9)]+[ss,ex],fasta=file_signature(fa),gtf=file_signature(gtf),prefix=str(prefix))
                log('Index complete: '+name)

    def extract_pilot(self,row,folder):
        one,two=folder/'pilot_R1.fastq',folder/'pilot_R2.fastq';n=0
        with open_fastq(row['r1']) as a,open_fastq(row['r2']) as b,open(one,'w') as x,open(two,'w') as y:
            for _ in range(self.c['pilot_pairs']):
                aa=[a.readline() for _ in range(4)];bb=[b.readline() for _ in range(4)]
                if not aa[0] and not bb[0]:break
                if not aa[-1] or not bb[-1] or read_id(aa[0])!=read_id(bb[0]):raise ValueError('Invalid/mismatched FASTQ pairs: '+row['sample'])
                x.writelines(aa);y.writelines(bb);n+=1
        if not n:raise ValueError('Empty pilot')
        return one,two,n

    def align_sort(self,row,index,strand,bam,tag,pilot=None):
        r1,r2=pilot if pilot else (row['r1'],row['r2'])
        cmd=[self.t['hisat2'],'-p',str(self.c['alignment_threads']),'--dta','--seed','42',*flags(strand,'hisat2'),'-x',index,'-1',r1,'-2',r2]
        if not pilot:cmd+=['--met-file',str(self.root/'logs'/f'{row["sample"]}.{tag}.metrics.tsv'),'--met','60']
        sorter=[self.t['samtools'],'sort','-@',str(self.c['sort_threads']),'-m',self.c.get('sort_memory','768M'),'-T',str(bam.parent/'sorttmp'),'-o',bam,'-']
        self.pipe(cmd,sorter,self.root/'logs'/f'{row["sample"]}.{tag}.hisat2.log',self.root/'logs'/f'{row["sample"]}.{tag}.sort.log')
        self.run([self.t['samtools'],'quickcheck',bam],self.root/'logs'/f'{row["sample"]}.samtools.log')

    def infer(self,row,scratch,index=None):
        sid=row['sample'];saved=self.root/'strand'/f'{sid}.json'
        if saved.exists():return json.loads(saved.read_text())
        if row.get('strandedness','auto')!='auto':
            result=dict(strand=row['strandedness'],status='user_supplied',informative_read1=0,forward_fraction=0,reverse_fraction=0)
            atomic_json(saved,result);return result
        one,two,n=self.extract_pilot(row,scratch);bam=scratch/'pilot.bam'
        self.align_sort(row,index or self.root/'index/extended','unstranded',bam,'pilot',pilot=(one,two))
        with open(self.root/'logs'/f'{sid}.strand.log','w') as err:
            process=subprocess.Popen([self.t['samtools'],'view','-F','2308','-f','64','-q','20',bam],stdout=subprocess.PIPE,stderr=err,text=True)
            counts=count_orientation(process.stdout,self.exons)
            process.stdout.close()
            if process.wait():raise RuntimeError('Strand SAM parsing failed: '+sid)
        result=decide_strand(counts['same'],counts['opposite'],min_reads=self.c['strand_min_reads'])
        result.update(pilot_pairs=n,ambiguous_exon_overlap=counts['ambiguous'],method='Unique MAPQ>=20 read1 overlaps reference exons; conflicting strands excluded.')
        atomic_json(saved,result)
        for path in (one,two,bam):self.safe_unlink(path)
        log(sid+': strand='+result['strand']+', status='+result['status']+', informative read1='+str(result['informative_read1']))
        return result

    def original_stats(self,row,strand):
        sid=row['sample'];path=self.root/'logs'/f'{sid}.original.flagstat.json'
        marker=self.root/'checkpoints'/f'{sid}.original.json'
        if self.checkpoint_valid(marker):
            json.loads(path.read_text());return
        log('Original-reference mapping: '+sid)
        align=[self.t['hisat2'],'-p',str(self.c['alignment_threads']),'--dta','--seed','42',*flags(strand,'hisat2'),
               '-x',self.root/'index/original','-1',row['r1'],'-2',row['r2'],'--met-file',self.root/'logs'/f'{sid}.original.metrics.tsv','--met','60']
        self.pipe(align,[self.t['samtools'],'flagstat','-@','2','-O','json','-'],self.root/'logs'/f'{sid}.original.hisat2.log',self.root/'logs'/f'{sid}.original.flagstat.log',stdout=path)
        json.loads(path.read_text())
        self.checkpoint(marker,[path,self.root/'logs'/f'{sid}.original.hisat2.log'],sample=sid)

    def first_pass(self,row):
        if self.stop.is_set():return
        sid=row['sample'];marker=self.root/'checkpoints'/f'{sid}.first.json';cram=self.root/'alignments'/f'{sid}.cram'
        assembly=self.root/'assembly'/f'{sid}.gtf'
        if self.checkpoint_valid(marker):
            self.run([self.t['samtools'],'quickcheck',cram],self.root/'logs'/f'{sid}.samtools.log');return
        self.space(int(row['r1_bytes'])+int(row['r2_bytes']))
        scratch=self.root/'scratch'/sid;scratch.mkdir(exist_ok=True)
        strand=self.infer(row,scratch)['strand'];bam=scratch/'full.bam'
        aligned=self.root/'checkpoints'/f'{sid}.aligned.json'
        if self.checkpoint_valid(aligned):
            self.run([self.t['samtools'],'view','-@',str(self.c['sort_threads']),'-b','-T',self.c['extended_fa'],'-o',bam,cram],self.root/'logs'/f'{sid}.decode.log')
            self.run([self.t['samtools'],'index','-@',str(self.c['sort_threads']),cram],self.root/'logs'/f'{sid}.cram.log')
        else:
            log('Extended-reference full mapping: '+sid)
            self.align_sort(row,self.root/'index/extended',strand,bam,'extended')
            self.run([self.t['samtools'],'flagstat','-O','json',bam],self.root/'logs'/f'{sid}.flagstat.log',stdout=self.root/'logs'/f'{sid}.extended.flagstat.json')
            self.run([self.t['samtools'],'view','-@',str(self.c['sort_threads']),'-C','--output-fmt-option','version=3.0','-T',self.c['extended_fa'],'-o',cram,bam],self.root/'logs'/f'{sid}.cram.log')
            self.run([self.t['samtools'],'quickcheck',cram],self.root/'logs'/f'{sid}.cram.log')
            self.run([self.t['samtools'],'index','-@',str(self.c['sort_threads']),cram],self.root/'logs'/f'{sid}.cram.log')
            self.checkpoint(aligned,[cram,Path(str(cram)+'.crai')],strand=strand)
        self.run([self.t['stringtie'],bam,'-p',str(self.c['stringtie_threads']),*flags(strand,'stringtie'),'-G',self.c['guide_gtf'],'-l',sid,'-o',assembly],self.root/'logs'/f'{sid}.assembly.log')
        nr_bam=self.root/'nonref_alignments'/f'{sid}.bam'
        self.run([self.t['samtools'],'view','-@','2','-b','-L',self.root/'annotation/nonref_regions.bed','-o',nr_bam,bam],self.root/'logs'/f'{sid}.nonref_bam.log')
        self.run([self.t['samtools'],'index',nr_bam],self.root/'logs'/f'{sid}.nonref_bam.log')
        self.safe_unlink(bam)
        if self.c['compare_original']:self.original_stats(row,strand)
        artifacts=[cram,Path(str(cram)+'.crai'),assembly,nr_bam,Path(str(nr_bam)+'.bai'),self.root/'strand'/f'{sid}.json']
        if self.c['compare_original']:artifacts+=[self.root/'logs'/f'{sid}.original.flagstat.json',self.root/'logs'/f'{sid}.original.hisat2.log']
        self.checkpoint(marker,artifacts,sample=sid,strand=strand)
        log('First pass complete: '+sid)

    def parallel(self,func,stage):
        atomic_json(self.root/'status.json',dict(state='RUNNING',stage=stage,samples=len(self.samples)))
        def reserved(row):
            if self.stop.is_set():return
            suffix='first' if stage=='mapping_and_assembly' else 'quant'
            if self.checkpoint_valid(self.root/'checkpoints'/f'{row["sample"]}.{suffix}.json'):
                return func(row)
            with self.reserve_disk(row):return func(row)
        pool=concurrent.futures.ThreadPoolExecutor(max_workers=self.c['concurrent_samples']);futures={pool.submit(reserved,row):row['sample'] for row in self.samples}
        done=0
        try:
            for future in concurrent.futures.as_completed(futures):
                future.result();done+=1;log(stage+': '+str(done)+'/'+str(len(self.samples))+' libraries complete')
        except BaseException:
            self.stop.set()
            for future in futures:future.cancel()
            pool.shutdown(wait=True,cancel_futures=True)
            raise
        else:pool.shutdown(wait=True)

    def merge(self):
        marker=self.root/'checkpoints/merged.json'
        if self.checkpoint_valid(marker):return
        # Rebuilding the common annotation invalidates every downstream quantification.
        for old in (self.root/'checkpoints').glob('*.quant.json'):self.safe_unlink(old)
        log('Merging observed RNA assemblies; no guide-only transcript insertion.')
        listing=self.root/'annotation/assembly_list.txt'
        listing.write_text(''.join(str(self.root/'assembly'/f'{r["sample"]}.gtf')+'\n' for r in self.samples))
        observed=self.root/'annotation/observed_merged.gtf'
        self.run([self.t['stringtie'],'--merge','-p',str(self.c['index_threads']),'-F','0','-T','0','-c','0','-l','BFRNA','-o',observed,listing],self.root/'logs/stringtie_merge.log')
        nrids=self.nonref_ids
        nr=self.root/'annotation/nonref.rna_merged.gtf';select_nonref_gtf(observed,nrids,nr)
        final=self.root/'annotation/final_annotation.gtf'
        with open(final,'w') as out:
            for source in (self.c['base_gtf'],nr):
                with open(source) as inp:
                    for line in inp:
                        if line.strip() and not line.startswith('#'):out.write(line.rstrip('\n')+'\n')
        gtf_genes(final,self.nonref_ids)
        self.run([self.t['gffread'],'-E',final,'-o','/dev/null'],self.root/'logs/final_annotation_validation.log')
        if nr.stat().st_size:
            self.run([self.t['gffcompare'],'-r',self.c['nonref_guide'],'-o',self.root/'annotation/nr_vs_augustus',nr],self.root/'logs/gffcompare.log')
            self.run([self.t['gffread'],nr,'-g',self.c['extended_fa'],'-w',self.root/'annotation/nonref.rna_transcripts.fa'],self.root/'logs/nonref_fasta.log')
        else:(self.root/'annotation/nonref.rna_transcripts.fa').write_text('')
        artifacts=[final,observed,nr,self.root/'annotation/nonref.rna_transcripts.fa',*list((self.root/'annotation').glob('nr_vs_augustus*.tmap'))]
        artifacts += [self.root/'assembly'/f'{r["sample"]}.gtf' for r in self.samples]
        self.checkpoint(marker,artifacts)

    def quantify(self,row):
        if self.stop.is_set():return
        sid=row['sample'];marker=self.root/'checkpoints'/f'{sid}.quant.json'
        quant=self.root/'quant'/sid;quant.mkdir(exist_ok=True);counts=self.root/'counts'/f'{sid}.tsv'
        if self.checkpoint_valid(marker):return
        self.space();scratch=self.root/'scratch'/sid;scratch.mkdir(exist_ok=True);bam=scratch/'quant.bam'
        strand=json.loads((self.root/'strand'/f'{sid}.json').read_text())['strand']
        log('Fixed-model quantification/counting: '+sid)
        self.run([self.t['samtools'],'view','-@',str(self.c['sort_threads']),'-b','-T',self.c['extended_fa'],'-o',bam,self.root/'alignments'/f'{sid}.cram'],self.root/'logs'/f'{sid}.quant_decode.log')
        self.run([self.t['stringtie'],bam,'-p',str(self.c['stringtie_threads']),*flags(strand,'stringtie'),'-e','-G',self.root/'annotation/final_annotation.gtf','-o',quant/'transcripts.gtf','-A',quant/'gene_abundance.tsv'],self.root/'logs'/f'{sid}.quant.log')
        self.run([self.t['featureCounts'],'-T',str(self.c['stringtie_threads']),'-p','--countReadPairs',*flags(strand,'featureCounts'),'-t','exon','-g','gene_id','-a',self.root/'annotation/final_annotation.gtf','-o',counts,bam],self.root/'logs'/f'{sid}.featureCounts.log')
        self.safe_unlink(bam)
        self.checkpoint(marker,[counts,quant/'gene_abundance.tsv',quant/'transcripts.gtf'],sample=sid)
        log('Quantification complete: '+sid)

    def run_all(self):
        self.preflight();self.build_indexes()
        nrbed=self.root/'annotation/nonref_regions.bed'
        with open(self.c['extended_fa']+'.fai') as inp,open(nrbed,'w') as out:
            for line in inp:
                c=line.split('\t')
                if c[0] in self.nonref_ids:out.write(f'{c[0]}\t0\t{c[1]}\n')
        log('Loading reference exon strand index.')
        self.exons=ExonIndex.from_gtf(self.c['base_gtf'])
        self.parallel(self.first_pass,'mapping_and_assembly')
        self.merge();self.parallel(self.quantify,'fixed_model_quantification')
        summary=collect(self.root,self.samples,self.c['core'],self.nonref_ids)
        atomic_json(self.root/'status.json',dict(state='COMPLETED',stage='results',summary=summary))
        (self.root/'COMPLETED').write_text(datetime.datetime.now().astimezone().isoformat()+'\n')
        log('COMPLETED: '+str(self.root/'results'))


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--config',required=True);parser.add_argument('--indexes-only',action='store_true');parser.add_argument('--wait-lock',action='store_true')
    args=parser.parse_args();config=json.loads(Path(args.config).read_text());runner=Pipeline(config)
    import fcntl
    lock=open(runner.root/'pipeline.lock','w')
    if args.wait_lock:log('Waiting for exclusive workspace lock (index builder may still be running).')
    try:fcntl.flock(lock,fcntl.LOCK_EX if args.wait_lock else fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise SystemExit('Another pipeline is already running in this workspace.')
    try:
        if args.indexes_only:runner.preflight();runner.build_indexes()
        else:runner.run_all()
    except BaseException as error:
        atomic_json(runner.root/'status.json',dict(state='FAILED',error=str(error),timestamp=datetime.datetime.now().astimezone().isoformat()))
        log('FAILED: '+str(error));raise


if __name__=='__main__':main()
