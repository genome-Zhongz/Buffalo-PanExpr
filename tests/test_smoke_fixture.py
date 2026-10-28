from pathlib import Path
import tempfile
import unittest
from make_smoke_fixture import make_fixture
from nonref_rna.config import load_config,resolve_reference
from nonref_rna.samples import load_samples
from nonref_rna.workflow import reference_contigs

class Synthetic(unittest.TestCase):
    def test_fixture_contains_unobserved_guide_and_six_paired_libraries(self):
        with tempfile.TemporaryDirectory() as d:
            p=make_fixture(Path(d)/'smoke')
            c=load_config(p/'config.yaml');refs=resolve_reference(c)
            self.assertEqual(len(load_samples(p/'samples.tsv')),6)
            self.assertEqual(reference_contigs(refs,c['reference']),['NR_1','NR_2','NR_unused'])
            self.assertTrue((p/'SYNTHETIC_DATA_ONLY.txt').exists())

if __name__=='__main__':unittest.main()
