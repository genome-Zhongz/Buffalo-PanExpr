import csv
from pathlib import Path
import tempfile
import unittest
import test_evidence
from nonref_rna.primary import export_primary_results
from nonref_rna.evidence import render_report

class Primary(unittest.TestCase):
    def test_primary_csv_and_fasta_are_P_candidates(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);test_evidence.Evidence().fixture(p);out=p/'export'
            de=out/'differential_expression';de.mkdir(parents=True)
            head='gene_id\tlogFC\tPValue\tFDR\tdirection\tis_nonreference\n'
            (de/'P05_FC1.nonreference_selected.tsv').write_text(head+'g1\t2\t0.01\t0.2\tUp\tTRUE\n')
            (de/'P05_FC1.all_selected.tsv').write_text(head+'g1\t2\t0.01\t0.2\tUp\tTRUE\n')
            result=export_primary_results(p,out,dict(expression=dict(tpm=100,min_samples=1)))
            self.assertEqual(result['nonreference_candidates'],1)
            self.assertEqual(result['rule'],'P < 0.05 and abs(log2FC) > 1')
            self.assertEqual(result['homology']['not_annotated'],1)
            self.assertTrue((out/'primary_Pvalue/homology_summary.md').is_file())
            with (out/'primary_Pvalue/nonreference_candidates_P05.csv').open(encoding='utf-8-sig') as f:rows=list(csv.DictReader(f))
            self.assertEqual(rows[0]['gene_id'],'g1')
            self.assertEqual(float(rows[0]['P']),0.01)
            self.assertNotIn('FDR',rows[0])
            self.assertIn('>g1 ',(out/'primary_Pvalue/nonreference_candidates_P05.longest.fa').read_text())
            h=out/'homology';h.mkdir()
            (h/'gene_homology.tsv').write_text('gene_id\tevidence_class\ng1\thigh_identity_broad_protein_coverage\nother\tno_detected_protein_match\n')
            result=export_primary_results(p,out,dict(expression=dict(tpm=100,min_samples=1)))
            self.assertEqual(result['homology']['with_matches'],1)
            self.assertEqual(result['homology']['candidates'],1)
            self.assertIn('1 个检出蛋白同源匹配',(out/'primary_Pvalue/homology_summary.md').read_text(encoding='utf-8'))
            render_report(out/'primary_Pvalue',{'primary_differential':result})
            html=(out/'primary_Pvalue/report/report.html').read_text(encoding='utf-8')
            self.assertIn('1 个检出蛋白同源匹配',html)
            self.assertNotIn('raw P',html)
            self.assertNotIn('exploratory',html)
            with (out/'primary_Pvalue/all_genes_P05.csv').open(encoding='utf-8-sig') as f:all_rows=list(csv.DictReader(f))
            self.assertEqual(float(all_rows[0]['P']),0.01)

if __name__=='__main__':unittest.main()
