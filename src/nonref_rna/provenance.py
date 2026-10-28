"""Check imported quantification, annotation and transcript reference identity."""
from collections import defaultdict
import json,math,re
from pathlib import Path
from .io import read_table,signature
from .engine.matrices import gtf_genes


def verify_annotation_and_quant(work,metadata,matrix,samples,nr_contigs,models,exons):
    work=Path(work);final=work/'annotation/final_annotation.gtf'
    if not final.is_file():raise ValueError('Import requires the final counting annotation GTF')
    actual=gtf_genes(final,nr_contigs)
    if set(actual)!=set(metadata):raise ValueError('Counting annotation and matrix gene IDs disagree')
    for g,m in metadata.items():
        a=actual[g]
        for key in ['contig','strand','start','end']:
            if str(m[key])!=str(a[key]):raise ValueError('Metadata/counting annotation disagreement: '+g+' '+key)
        if (str(m['is_nonreference']).lower() in ('1','true'))!=a['is_nonreference']:raise ValueError('Nonreference flag/counting annotation disagreement: '+g)
    counted_models={};counted_exons=defaultdict(list)
    with final.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'):continue
            c=line.rstrip().split('\t');a=dict(re.findall(r'(\w+) "([^"]*)"',c[8]));tid=a.get('transcript_id')
            if tid not in models:continue
            if c[2]=='transcript':counted_models[tid]=dict(gene_id=a['gene_id'],contig=c[0],strand=c[6],start=int(c[3]),end=int(c[4]))
            if c[2]=='exon':counted_exons[tid].append((int(c[3]),int(c[4])))
    for tid,model in models.items():
        expected_model={k:model[k] for k in ['gene_id','contig','strand','start','end']}
        if counted_models.get(tid)!=expected_model or set(counted_exons[tid])!=set(exons[tid]):
            raise ValueError('Exported transcript/counting annotation structure disagreement: '+tid)
    expected={r['gene_id']:r for r in matrix}
    for sample in samples:
        sid=sample['sample'];path=work/'quant'/sid/'gene_abundance.tsv'
        if not path.is_file():raise ValueError('Import requires source StringTie quantification: '+str(path))
        values=defaultdict(float)
        for r in read_table(path):
            g=r.get('Gene ID')
            if not g or g in ('-','.'):continue
            if g not in metadata:raise ValueError('Source quantification contains foreign gene ID: '+g)
            m=metadata[g]
            if r.get('Reference')!=m['contig']:raise ValueError('Quantification/annotation contig disagreement: '+g)
            if r.get('Strand')!=m['strand']:raise ValueError('Quantification/annotation strand disagreement: '+g)
            if not int(m['start'])<=int(r['Start'])<=int(r['End'])<=int(m['end']):raise ValueError('Quantification outside annotated locus: '+g)
            values[g]+=float(r['TPM'])
        for g in metadata:
            v=float(expected[g][sid])
            if not math.isclose(values[g],v,rel_tol=6e-8,abs_tol=1e-7):raise ValueError('TPM matrix/source quantification disagreement: '+g+' '+sid)


def reference_from_manifest(work):
    work=Path(work);config=None
    for name in ['run_config.snapshot.json','config.json']:
        p=work/name
        if p.is_file():
            c=json.loads(p.read_text(encoding='utf-8'))
            if c.get('extended_fa'):config=c;break
    if config is None:raise ValueError('Import requires source run config with extended_fa')
    path=Path(config['extended_fa'])
    if not path.is_file():raise ValueError('Source extended reference no longer available: '+str(path))
    expected=[]
    for name in ['run_manifest.json','analysis_identity.json']:
        p=work/name
        if p.is_file():
            identity=json.loads(p.read_text(encoding='utf-8')).get('identity',{})
            expected+=identity.get('inputs',[])+identity.get('references',[])
    current=signature(path)
    if current not in expected:raise ValueError('Extended reference is not verified by the source input manifest')
    return path


class ReferenceSequence:
    def __init__(self,path):
        self.handle=open(path,'rb');self.index={};idx=Path(str(path)+'.fai')
        if idx.is_file():
            for r in idx.read_text().splitlines():
                c=r.split('\t');self.index[c[0]]=tuple(map(int,c[1:5]))
        else:
            # Small fixtures and references without faidx: collect offsets only.
            ident=None
            while True:
                pos=self.handle.tell();line=self.handle.readline()
                if not line:break
                if line.startswith(b'>'):
                    ident=line[1:].split()[0].decode();self.index[ident]=[0,self.handle.tell(),0,0]
                elif line.strip():
                    if ident is None:raise ValueError('Invalid source FASTA')
                    record=self.index[ident];record[0]+=len(line.rstrip(b'\r\n'))
                    if not record[2]:record[2]=len(line.rstrip(b'\r\n'));record[3]=len(line)

    def fetch(self,chrom,start,end):
        if chrom not in self.index:raise ValueError('Unknown source reference contig '+chrom)
        length,offset,bases,width=self.index[chrom]
        if not 1<=start<=end<=length or not bases:raise ValueError('Sequence outside reference')
        a=start-1;b=end-1
        first=offset+(a//bases)*width+a%bases
        last=offset+(b//bases)*width+b%bases+1
        self.handle.seek(first)
        return self.handle.read(last-first).replace(b'\n',b'').replace(b'\r',b'').decode('ascii')

    def close(self):self.handle.close()


def verify_transcript_sequences(work,models,exons,sequences):
    path=reference_from_manifest(work);reader=ReferenceSequence(path)
    complement=str.maketrans('ACGTRYSWKMBDHVN','TGCAYRSWMKVHDBN')
    try:
        for tid,m in models.items():
            intervals=[]
            for a,b in sorted(set(exons[tid])):
                if intervals and a<=intervals[-1][1]+1:intervals[-1][1]=max(intervals[-1][1],b)
                else:intervals.append([a,b])
            seq=''.join(reader.fetch(m['contig'],a,b) for a,b in intervals).upper()
            if m['strand']=='-':seq=seq.translate(complement)[::-1]
            if seq!=sequences[tid].upper():raise ValueError('Transcript/reference sequence disagreement: '+tid)
    finally:reader.close()
    return str(path)
