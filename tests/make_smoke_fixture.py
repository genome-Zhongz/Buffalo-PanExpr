"""Generate synthetic DNA/RNA reads solely for external-tool integration tests."""
import argparse,csv,json,random
from pathlib import Path
import yaml


def make_fixture(folder,tools=None):
    p=Path(folder).resolve()
    if p.exists():raise ValueError('Use a new fixture directory')
    p.mkdir(parents=True);core=p/'bundle'
    for sub in ['rnaseq_reference','reference','reads']:(core/sub).mkdir(parents=True)
    rng=random.Random(71411)
    seqs={k:''.join(rng.choices('ACGT',k=n)) for k,n in [('chr1',10000),('NR_1',6000),('NR_2',5000),('NR_unused',3000)]}
    def fa(path,names):
        with path.open('w') as f:
            for k in names:f.write('>'+k+'\n'+seqs[k]+'\n')
    fa(core/'rnaseq_reference/extended.fa',seqs)
    fa(core/'reference/original.fa',['chr1'])
    genes=[('refA','chr1',1000,2200),('refB','chr1',4000,5200),('pred1','NR_1',800,2000),('pred2','NR_2',1000,2200),('unused','NR_unused',500,1700)]
    def gtf(rows):
        return ''.join(f'{c}\tsynthetic\t{ft}\t{a}\t{b}\t.\t+\t.\tgene_id "{g}"; transcript_id "{g}.1";\n' for g,c,a,b in rows for ft in ['transcript','exon'])
    (core/'reference/reference.normalized.gtf').write_text(gtf(genes[:2]))
    (core/'rnaseq_reference/nonref.guide.gtf').write_text(gtf(genes[2:]))
    (core/'rnaseq_reference/guide.gtf').write_text(gtf(genes))
    (core/'run_config.json').write_text(json.dumps(dict(graph_reference=str(core/'reference/original.fa'))))
    complement=str.maketrans('ACGT','TGCA');rows=[]
    for i in range(6):
        sid=f'synthetic{i+1}';r1=core/'reads'/(sid+'_R1.fastq');r2=core/'reads'/(sid+'_R2.fastq')
        with r1.open('w') as f1,r2.open('w') as f2:
            counter=0
            amounts=[800+rng.randrange(70),700+rng.randrange(70),(300 if i<3 else 650)+rng.randrange(50),400+rng.randrange(50)]
            for (_,chrom,a,b),number in zip(genes,amounts):
                template=seqs[chrom][a-1:b]
                for _ in range(number):
                    start=rng.randrange(len(template)-300);fragment=template[start:start+300]
                    one=fragment[:100];two=fragment[-100:].translate(complement)[::-1]
                    counter+=1;name=f'{sid}_{counter}'
                    f1.write(f'@{name}/1\n{one}\n+\n'+('I'*100)+'\n')
                    f2.write(f'@{name}/2\n{two}\n+\n'+('I'*100)+'\n')
        rows.append(dict(sample=sid,r1=str(r1),r2=str(r2),tissue='Synthetic',breed='Synthetic',group='Control' if i<3 else 'Treatment',individual=sid,strandedness='unstranded'))
    with (p/'samples.tsv').open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
    cfg=dict(reference=dict(bundle=str(core)),reads=dict(cleaned=False),resources=dict(threads=8,jobs=1,minimum_free_gib=0),
             expression=dict(tpm=1,min_samples=2),analysis=dict(plots=True,differential=True,control='Control',treatment='Treatment',cpm=1,min_samples=2),
             homology=dict(enabled=False),tools=tools or {})
    (p/'config.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    (p/'SYNTHETIC_DATA_ONLY.txt').write_text('Random synthetic sequences and reads. No biological findings may be inferred. Expected: two observed nonreference loci; unused guide must not become expressed candidate.\n')
    return p


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('out');parser.add_argument('--tools-json')
    a=parser.parse_args();tools=json.loads(Path(a.tools_json).read_text()) if a.tools_json else None
    print(make_fixture(a.out,tools))
