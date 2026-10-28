"""Public CLI. Full alignment runs require Linux; summaries are portable."""
import argparse
import copy
import json
from pathlib import Path
import sys
import subprocess
import yaml
from . import __version__
from .config import DEFAULTS,load_config


def main(argv=None):
    p=argparse.ArgumentParser(prog='nonref-rna',description='RNA-supported nonreference gene candidates from paired-end RNA-seq and a prepared reference bundle.')
    p.add_argument('--version',action='version',version=__version__)
    commands=p.add_subparsers(dest='command',required=True)
    init=commands.add_parser('init',help='Create configuration for a prepared reference bundle')
    init.add_argument('--reference-bundle',required=True);init.add_argument('--out',required=True)
    for name in ['validate','run']:
        q=commands.add_parser(name);q.add_argument('--config',required=True);q.add_argument('--samples',required=True);q.add_argument('--out',default='nonref_results')
        if name=='run':q.add_argument('--dry-run',action='store_true')
    q=commands.add_parser('summarize',help='Import a compatible completed RNA workspace and export candidate evidence')
    q.add_argument('--work',required=True);q.add_argument('--out',required=True);q.add_argument('--config')
    q=commands.add_parser('status');q.add_argument('--out',required=True)
    a=p.parse_args(argv)
    try:
        if a.command=='init':
            target=Path(a.out).resolve()
            if target.exists():raise ValueError('Refusing to overwrite existing config: '+str(target))
            core=Path(a.reference_bundle).expanduser().resolve()
            if not core.is_dir():raise ValueError('Reference bundle directory does not exist')
            c=copy.deepcopy(DEFAULTS);c['reference']['bundle']=str(core)
            db=core.parent/'combined_db.dmnd'
            if db.is_file():c['homology']['database']=str(db)
            target.parent.mkdir(parents=True,exist_ok=True)
            target.write_text(yaml.safe_dump(c,sort_keys=False),encoding='utf-8')
            sample=target.parent/'samples.example.tsv'
            if not sample.exists():sample.write_text('sample\tr1\tr2\ttissue\tbreed\tgroup\tindividual\tstrandedness\nsample1\t/path/to/R1.fastq.gz\t/path/to/R2.fastq.gz\tMammary\tBuffalo\tEarly\tindividual1\tauto\n')
            print('Created '+str(target));return 0
        if a.command=='status':
            source=Path(a.out)/'status.json'
            if not source.is_file():raise ValueError('No status.json in '+a.out)
            print(source.read_text(encoding='utf-8'));return 0
        from .workflow import validate,run,postprocess
        if a.command=='summarize':
            c=load_config(a.config) if a.config else copy.deepcopy(DEFAULTS)
            if Path(a.work).resolve()==Path(a.out).resolve():raise ValueError('For imported results use a NEW export directory, not the source workspace')
            if (Path(a.out)/'report/summary.json').exists():raise ValueError('Export already exists; use a new output directory')
            result=postprocess(c,a.work,a.out)
        else:
            c=load_config(a.config)
            if a.command=='validate':
                rows,engine,nr=validate(c,a.samples,a.out)
                result=dict(valid=True,libraries=len(rows),nonreference_contigs=len(nr),tools=engine['tools'],resources={k:engine[k] for k in ['index_threads','alignment_threads','sort_threads','concurrent_samples']})
            else:result=run(c,a.samples,a.out,a.dry_run)
        print(json.dumps(result,indent=2,ensure_ascii=False));return 0
    except (OSError,ValueError,RuntimeError,KeyError,subprocess.CalledProcessError) as e:
        print('ERROR: '+str(e),file=sys.stderr);return 2
