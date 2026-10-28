#!/usr/bin/env python3
import bisect
import re


def merged(intervals):
    result=[]
    for a,b in sorted(intervals):
        if result and a<=result[-1][1]:result[-1][1]=max(result[-1][1],b)
        else:result.append([a,b])
    return result


def cigar_blocks(position,cigar):
    tokens=re.findall(r'(\d+)([MIDNSHP=X])',cigar)
    if ''.join(n+op for n,op in tokens)!=cigar:raise ValueError('Invalid CIGAR: '+cigar)
    pos=position-1; blocks=[]
    for n,op in tokens:
        size=int(n)
        if op in 'M=X':blocks.append((pos,pos+size));pos+=size
        elif op in 'DN':pos+=size
    return blocks


class ExonIndex:
    def __init__(self,intervals):
        self.data={}
        for key,values in intervals.items():
            data=merged(values)
            self.data[key]=([a for a,b in data],[b for a,b in data])

    @classmethod
    def from_gtf(cls,path):
        from collections import defaultdict
        intervals=defaultdict(list)
        with open(path) as inp:
            for line in inp:
                if line.startswith('#') or not line.strip():continue
                c=line.rstrip('\n').split('\t')
                if c[2]=='exon' and c[6] in ('+','-'):
                    intervals[(c[0],c[6])].append((int(c[3])-1,int(c[4])))
        return cls(intervals)

    def strands(self,chrom,blocks):
        found=set()
        for strand in ('+','-'):
            data=self.data.get((chrom,strand))
            if not data:continue
            starts,ends=data
            for a,b in blocks:
                i=bisect.bisect_left(starts,b)-1
                if i>=0 and ends[i]>a:found.add(strand);break
        return found


def count_orientation(lines,index):
    counts=dict(same=0,opposite=0,ambiguous=0)
    for line in lines:
        if line.startswith('@'):continue
        c=line.rstrip('\n').split('\t');flag=int(c[1])
        if flag&2308 or not flag&64 or int(c[4])<20 or 'NH:i:1' not in c[11:]:continue
        strands=index.strands(c[2],cigar_blocks(int(c[3]),c[5]))
        if len(strands)!=1:counts['ambiguous']+=1;continue
        gene_reverse=next(iter(strands))=='-'; read_reverse=bool(flag&16)
        counts['same' if read_reverse==gene_reverse else 'opposite']+=1
    return counts


def decide_strand(same,opposite,min_reads=1000,threshold=.8):
    total=same+opposite
    result=dict(informative_read1=total,same_gene_strand=same,opposite_gene_strand=opposite,
                forward_fraction=same/total if total else 0,reverse_fraction=opposite/total if total else 0)
    if total<min_reads:strand,status='unstranded','insufficient_informative_reads'
    elif opposite/total>=threshold:strand,status='rf','high_confidence'
    elif same/total>=threshold:strand,status='fr','high_confidence'
    else:strand,status='unstranded','mixed_or_unstranded'
    result.update(strand=strand,status=status)
    return result
