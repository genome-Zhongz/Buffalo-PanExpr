import unittest
from nonref_rna.homology import classify_hit, summarize_homology

class Homology(unittest.TestCase):
    def test_candidate_summary_counts_unique_genes_and_unannotated_separately(self):
        classes=(['high_identity_broad_protein_coverage']*25+
                 ['high_identity_partial_protein_match']*12+
                 ['limited_or_divergent_protein_match']*5+
                 ['no_detected_protein_match']*6)
        rows=[dict(gene_id=str(i),evidence_class=c) for i,c in enumerate(classes)]
        rows.append(dict(gene_id='outside',evidence_class=classes[0]))
        result=summarize_homology(rows,[str(i) for i in range(48)])
        self.assertEqual(result['with_matches'],42)
        self.assertEqual(result['high_identity_broad_protein_coverage'],25)
        self.assertEqual(result['without_matches'],6)
        self.assertIn('48 个候选中，42 个检出蛋白同源匹配',result['text_zh'])
        missing=summarize_homology([],['g1'])
        self.assertEqual(missing['not_annotated'],1)
        self.assertEqual(missing['without_matches'],0)
        self.assertIn('未注释',missing['text_zh'])

    def test_subject_coverage_and_short_hit_not_overcalled(self):
        self.assertEqual(classify_hit({'pident':99,'scovhsp':100,'length':120,'evalue':1e-30}),'high_identity_broad_protein_coverage')
        self.assertEqual(classify_hit({'pident':99,'scovhsp':20,'length':120,'evalue':1e-30}),'high_identity_partial_protein_match')
        self.assertEqual(classify_hit({'pident':100,'scovhsp':100,'length':20,'evalue':1e-6}),'limited_or_divergent_protein_match')

if __name__=='__main__':unittest.main()
