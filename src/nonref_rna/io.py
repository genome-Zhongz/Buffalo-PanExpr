import csv
import gzip
import hashlib
import json
import os
from pathlib import Path


def open_fastq(path):
    return gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)


def fasta(path):
    name=None;seq=[]
    with open(path) as f:
        for line in f:
            if line.startswith('>'):
                if name is not None:yield name,''.join(seq)
                name=line[1:].strip();seq=[]
            elif line.strip():
                if name is None:raise ValueError('Sequence before FASTA header: '+str(path))
                seq.append(line.strip())
    if name is not None:yield name,''.join(seq)


def write_fasta(handle,header,seq):
    handle.write('>'+header+'\n')
    for i in range(0,len(seq),80):handle.write(seq[i:i+80]+'\n')


def read_table(path):
    csv.field_size_limit(100000000)
    with open(path,encoding='utf-8-sig') as f:return list(csv.DictReader(f,delimiter='\t'))


def write_table(path,fields,rows):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore',lineterminator='\n')
        w.writeheader();w.writerows(rows)


def atomic_json(path,data):
    path=Path(path);path.parent.mkdir(parents=True,exist_ok=True)
    tmp=path.with_suffix(path.suffix+'.tmp')
    tmp.write_text(json.dumps(data,indent=2,ensure_ascii=False)+'\n',encoding='utf-8');os.replace(tmp,path)


def signature(path):
    p=Path(path).resolve();s=p.stat()
    return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns)


def digest(path):
    h=hashlib.sha256()
    with open(path,'rb') as f:
        for block in iter(lambda:f.read(2**20),b''):h.update(block)
    return h.hexdigest()
