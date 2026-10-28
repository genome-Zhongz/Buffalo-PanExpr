import contextlib
import io
from pathlib import Path
import tempfile
import unittest
from nonref_rna.cli import main

class CLI(unittest.TestCase):
    def test_init_generates_portable_config_and_does_not_overwrite(self):
        with tempfile.TemporaryDirectory() as d:
            p=Path(d);bundle=p/'bundle';bundle.mkdir()
            with contextlib.redirect_stdout(io.StringIO()):
                self.assertEqual(main(['init','--reference-bundle',str(bundle),'--out',str(p/'config.yaml')]),0)
                self.assertNotEqual(main(['init','--reference-bundle',str(bundle),'--out',str(p/'config.yaml')]),0)
            text=(p/'config.yaml').read_text(encoding='utf-8')
            self.assertIn('bundle:',text)
            self.assertNotIn('sshpass',text)

if __name__=='__main__':unittest.main()
