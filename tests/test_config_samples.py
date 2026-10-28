import gzip
from pathlib import Path
import tempfile
import unittest

from nonref_rna.config import build_resources, load_config
from nonref_rna.samples import load_samples


class Inputs(unittest.TestCase):
    def test_resource_budget(self):
        r=build_resources(64,3)
        self.assertLessEqual(r['concurrent_samples']*(r['alignment_threads']+r['sort_threads']+1),64)
        self.assertEqual(r['index_threads'],64)
        with self.assertRaises(ValueError):build_resources(2,3)

    def reads(self,p):
        (p/'one.fastq').write_text('@read/1\nACGT\n+\nIIII\n')
        with gzip.open(p/'two.fastq.gz','wt') as f:f.write('@read/2\nACGT\n+\nIIII\n')

    def test_relative_plain_and_gzip_reads(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.reads(p)
            (p/'samples.tsv').write_text('sample\tr1\tr2\ns1\tone.fastq\ttwo.fastq.gz\n')
            rows=load_samples(p/'samples.tsv')
            self.assertEqual(rows[0]['sample'],'s1')
            self.assertEqual(rows[0]['r1'],str((p/'one.fastq').resolve()))
            self.assertEqual(rows[0]['strandedness'],'auto')

    def test_duplicate_and_unsafe_sample_id_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.reads(p)
            for body in ['s1\tone.fastq\ttwo.fastq.gz\ns1\tone.fastq\ttwo.fastq.gz\n','../x\tone.fastq\ttwo.fastq.gz\n']:
                (p/'samples.tsv').write_text('sample\tr1\tr2\n'+body)
                with self.assertRaises(ValueError):load_samples(p/'samples.tsv')

    def test_invalid_threshold_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)/'config.yaml';p.write_text('reference:\n  bundle: .\nexpression:\n  min_samples: 0\n')
            with self.assertRaises(ValueError):load_config(p)

if __name__=='__main__':unittest.main()
