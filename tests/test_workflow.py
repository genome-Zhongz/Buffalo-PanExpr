from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch
import types,json,sys
from nonref_rna.config import DEFAULTS
from nonref_rna.workflow import run
from nonref_rna.workflow import check_identity, reference_contigs, prepare_engine_config, validate_analysis
from nonref_rna.engine.pipeline import Pipeline, splice_index_options


class Workflow(unittest.TestCase):
    def test_incompatible_resume_preserves_completed_status(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);ref=p/'reference';ref.write_text('FASTA')
            (p/'COMPLETED').write_text('done')
            (p/'status.json').write_text('{"state":"COMPLETED"}')
            (p/'run_manifest.json').write_text('{"fingerprint":"different"}')
            engine={k:str(ref) for k in ['extended_fa','guide_gtf','base_gtf','nonref_guide','base_fa']};engine['tools']={}
            rows=[dict(sample='s1',r1=str(ref),r2=str(ref))]
            fake=types.SimpleNamespace(LOCK_EX=1,LOCK_NB=2,flock=lambda *a:None)
            with patch('nonref_rna.workflow.validate',return_value=(rows,engine,['NR_1'])),patch('nonref_rna.workflow.sys.platform','linux'),patch.dict(sys.modules,{'fcntl':fake}):
                with self.assertRaises(ValueError):run(DEFAULTS,'unused',p)
            self.assertTrue((p/'COMPLETED').exists())
            self.assertEqual(json.loads((p/'status.json').read_text())['state'],'COMPLETED')

    def test_repeated_individual_rejected_for_unpaired_model(self):
        rows=[{'sample':str(i),'group':g,'individual':animal} for i,(g,animal) in enumerate([('A','x'),('A','y'),('B','x'),('B','z')])]
        with self.assertRaises(ValueError):validate_analysis(rows,{'differential':True,'control':'A','treatment':'B','min_samples':2})
        rows[2]['individual']='w'
        validate_analysis(rows,{'differential':True,'control':'A','treatment':'B','min_samples':2})

    def test_single_exon_guide_uses_linear_index_without_empty_splice_options(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);ss=p/'sites';ex=p/'exons';ss.write_text('');ex.write_text('chr1\t0\t999\n')
            self.assertEqual(splice_index_options(ss,ex),[])
            ss.write_text('chr1\t100\t200\t+\n')
            self.assertEqual(splice_index_options(ss,ex),['--ss',str(ss),'--exon',str(ex)])

    def test_identity_allows_resources_but_not_input_change(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);source=p/'reads';source.write_text('original')
            state={'inputs':[str(source)],'biology':{'cleaned':True},'resources':{'threads':8}}
            check_identity(p,state)
            state['resources']['threads']=16;check_identity(p,state)
            source.write_text('changed input')
            with self.assertRaises(ValueError):check_identity(p,state)

    def test_nonreference_selection_and_guide_validation(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d)
            (p/'extended.fa').write_text('>chr1\nACGT\n>NR_1\nACGT\n')
            (p/'base.gtf').write_text('chr1\tx\texon\t1\t4\t.\t+\t.\tgene_id "ref"; transcript_id "ref.1";\n')
            (p/'nr.gtf').write_text('NR_1\tx\texon\t1\t4\t.\t+\t.\tgene_id "pred"; transcript_id "pred.1";\n')
            cfg={'extended_fa':str(p/'extended.fa'),'base_gtf':str(p/'base.gtf'),'guide_gtf':str(p/'nr.gtf'),'nonref_guide':str(p/'nr.gtf')}
            result=reference_contigs(cfg,{'nonref_prefix':'NR_'})
            self.assertEqual(result,['NR_1'])
            (p/'nr.gtf').write_text('chrX\tx\texon\t1\t4\t.\t+\t.\tgene_id "bad"; transcript_id "bad.1";\n')
            with self.assertRaises(ValueError):reference_contigs(cfg,{'nonref_prefix':'NR_'})

    def test_checkpoint_detects_edited_artifact(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);(p/'samples.tsv').write_text('sample\ns1\n')
            runner=Pipeline({'work':d,'tools':{}})
            a=p/'artifact';a.write_text('ok');marker=p/'checkpoints/test.json'
            runner.checkpoint(marker,[a]);self.assertTrue(runner.checkpoint_valid(marker))
            a.write_text('edited');self.assertFalse(runner.checkpoint_valid(marker))

if __name__=='__main__':unittest.main()
