import csv
from pathlib import Path
import re
from .io import open_fastq,write_table


def read_id(header):
    token=header.strip().split()[0]
    return token[:-2] if token.endswith(('/1','/2')) else token


def first_record(path):
    with open_fastq(path) as f:a=[f.readline().rstrip('\n\r') for _ in range(4)]
    if not a[0].startswith('@') or not a[2].startswith('+') or not a[1] or len(a[1])!=len(a[3]):
        raise ValueError('Invalid first FASTQ record: '+str(path))
    return a


def load_samples(path):
    path=Path(path).resolve()
    with path.open(encoding='utf-8-sig') as f:
        reader=csv.DictReader(f,delimiter='\t')
        if not {'sample','r1','r2'}<=set(reader.fieldnames or []):raise ValueError('Samples require sample, r1, r2 columns')
        rows=list(reader)
    if not rows:raise ValueError('Samples table is empty')
    seen=set();used=set()
    for row in rows:
        sid=row['sample']
        if not re.fullmatch('[A-Za-z0-9][A-Za-z0-9_.-]*',sid or '') or sid in seen:raise ValueError('Duplicate/unsafe sample ID: '+str(sid))
        seen.add(sid)
        for k in ['tissue','breed','group','individual']:row[k]=row.get(k) or 'unspecified'
        row['strandedness']=row.get('strandedness') or 'auto'
        if row['strandedness'] not in ['auto','rf','fr','unstranded']:raise ValueError('Invalid strandedness for '+sid)
        for k in ['r1','r2']:
            p=Path(row[k]).expanduser();p=(path.parent/p if not p.is_absolute() else p).resolve()
            if not p.is_file() or not p.stat().st_size:raise ValueError('Missing/empty reads: '+str(p))
            if p in used:raise ValueError('FASTQ file reused across mates/samples: '+str(p))
            used.add(p);row[k]=str(p);row[k+'_bytes']=p.stat().st_size
        a,b=first_record(row['r1']),first_record(row['r2'])
        if read_id(a[0])!=read_id(b[0]):raise ValueError('Read mate IDs disagree: '+sid)
    return rows


def save_samples(path,rows):
    fields=['sample','r1','r2','tissue','breed','group','individual','strandedness','r1_bytes','r2_bytes']
    write_table(path,fields,rows)
