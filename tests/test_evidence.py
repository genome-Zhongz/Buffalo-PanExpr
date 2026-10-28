import csv,json
from pathlib import Path
import tempfile
import unittest
from nonref_rna.evidence import summarize
from nonref_rna.io import signature


class Evidence(unittest.TestCase):
    def fixture(self,p):
        for d in ['results','annotation']:(p/d).mkdir()
        (p/'samples.tsv').write_text('sample\ttissue\tbreed\ns1\tBrain\tA\ns2\tBrain\tA\ns3\tHeart\tB\n')
        (p/'results/gene_TPM.tsv').write_text('gene_id\ts1\ts2\ts3\nref\t4\t4\t0\ng1\t1\t0\t0\ng2\t3\t3\t0\ng3\t0\t0\t0\n')
        (p/'results/gene_metadata.tsv').write_text('gene_id\tcontig\tstrand\tstart\tend\tgene_name\tis_nonreference\nref\tchr1\t+\t1\t5\tx\t0\ng1\tNR_1\t+\t1\t6\t\t1\ng2\tNR_2\t+\t1\t9\t\t1\ng3\tNR_3\t+\t1\t3\t\t1\n')
        lines=[]
        for g,t,c,a,b in [('g1','g1.1','NR_1',1,3),('g1','g1.2','NR_1',1,6),('g2','g2.1','NR_2',1,9),('g3','g3.1','NR_3',1,3)]:
            for ft in ['transcript','exon']:lines.append(f'{c}\tStringTie\t{ft}\t{a}\t{b}\t.\t+\t.\tgene_id "{g}"; transcript_id "{t}";\n')
        (p/'annotation/nonref.rna_merged.gtf').write_text(''.join(lines))
        (p/'annotation/nonref.rna_transcripts.fa').write_text('>g1.1\nATG\n>g1.2\nATGAAA\n>g2.1\nATGAAATAA\n>g3.1\nATG\n')
        (p/'annotation/final_annotation.gtf').write_text('chr1\tx\texon\t1\t5\t.\t+\t.\tgene_id "ref"; transcript_id "ref.1";\n'+''.join(lines))
        reference=p/'extended.fa';reference.write_text('>chr1\nAAAAA\n>NR_1\nATGAAA\n>NR_2\nATGAAATAA\n>NR_3\nATG\n')
        (p/'run_config.snapshot.json').write_text(json.dumps(dict(extended_fa=str(reference))))
        (p/'run_manifest.json').write_text(json.dumps(dict(identity=dict(inputs=[signature(reference)]))))
        values={'s1':[4,1,3,0],'s2':[4,0,3,0],'s3':[0,0,0,0]}
        for sample,v in values.items():
            q=p/'quant'/sample;q.mkdir(parents=True)
            data='Gene ID\tReference\tStrand\tStart\tEnd\tTPM\n'
            for gene,contig,end,tpm in zip(['ref','g1','g2','g3'],['chr1','NR_1','NR_2','NR_3'],[5,6,9,3],v):data+=f'{gene}\t{contig}\t+\t1\t{end}\t{tpm}\n'
            (q/'gene_abundance.tsv').write_text(data)

    def test_foreign_metadata_and_quantification_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            f=p/'results/gene_metadata.tsv';f.write_text(f.read_text().replace('g2\tNR_2\t+','g2\tNR_2\t-'))
            with self.assertRaises(ValueError):summarize(p,p/'export',dict(expression=dict(tpm=1,min_samples=1)))
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            f=p/'results/gene_TPM.tsv';f.write_text(f.read_text().replace('g2\t3\t3\t0','g2\t4\t3\t0'))
            with self.assertRaises(ValueError):summarize(p,p/'export',dict(expression=dict(tpm=1,min_samples=1)))

    def test_wrong_sequence_same_length_rejected(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            f=p/'annotation/nonref.rna_transcripts.fa';f.write_text(f.read_text().replace('ATGAAATAA','ATGCCCTAA'))
            with self.assertRaises(ValueError):summarize(p,p/'export',dict(expression=dict(tpm=1,min_samples=1)))

    def test_exported_model_must_match_counting_annotation_strand(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            f=p/'annotation/nonref.rna_merged.gtf'
            f.write_text(f.read_text().replace('NR_2\tStringTie\ttranscript\t1\t9\t.\t+','NR_2\tStringTie\ttranscript\t1\t9\t.\t-').replace('NR_2\tStringTie\texon\t1\t9\t.\t+','NR_2\tStringTie\texon\t1\t9\t.\t-'))
            f=p/'annotation/nonref.rna_transcripts.fa';f.write_text(f.read_text().replace('ATGAAATAA','TTATTTCAT'))
            with self.assertRaises(ValueError):summarize(p,p/'export',dict(expression=dict(tpm=1,min_samples=1)))

    def test_candidate_dedup_and_all_isoforms(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            report=summarize(p,p/'export',dict(expression=dict(tpm=1,min_samples=1)))
            self.assertEqual(report['candidate_genes'],2)
            self.assertEqual(report['candidate_transcripts'],3)
            self.assertEqual(report['all_genes_tissue_mean_gt1_union'],2)
            self.assertEqual(report['nonref_genes_tissue_mean_gt1_union'],1)
            fa=(p/'export/candidates/longest_transcript_per_gene.fa').read_text()
            self.assertIn('>g1 transcript_id=g1.2',fa)
            self.assertNotIn('>ref',fa)
            self.assertEqual(report['homology_status'],'not_run')

    def test_empty_set_has_valid_files_and_summary(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            report=summarize(p,p/'export',dict(expression=dict(tpm=100,min_samples=1)))
            self.assertEqual(report['candidate_genes'],0)
            self.assertEqual((p/'export/candidates/longest_transcript_per_gene.fa').read_text(),'')

    def test_explicit_differential_gene_list_not_filtered_again_by_TPM(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            report=summarize(p,p/'export',dict(expression=dict(tpm=100,min_samples=1)),gene_ids=['g1'],copy_matrices=False)
            self.assertEqual(report['candidate_genes'],1)
            self.assertFalse((p/'export/matrices').exists())

    def test_missing_transcript_is_error(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);self.fixture(p)
            (p/'annotation/nonref.rna_transcripts.fa').write_text('>g1.1\nATG\n')
            with self.assertRaises(ValueError):summarize(p,p/'export',dict(expression=dict(tpm=1,min_samples=1)))

if __name__=='__main__':unittest.main()
