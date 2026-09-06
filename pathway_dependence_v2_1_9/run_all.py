"""Verify the frozen archive, reconstruct in isolation, compare, then reverify."""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tempfile
sys.dont_write_bytecode=True
from compare_reconstruction import compare_file

ROOT=Path(__file__).resolve().parent
MANUSCRIPT_DERIVED=(
    '37_REAL_CFD_COMMON_PARTITION_PROFILES.csv',
    '37_REAL_CFD_COMMON_PARTITION_PROFILES.json',
    '42_REAL_CFD_COMMON_PARTITION_SUMMARY.csv',
    '42_REAL_CFD_COMMON_PARTITION_SUMMARY.json',
    'weight_stability_results.json',
)

AUXILIARY_DERIVED=('auxiliary/robust_union/28_ROBUST_UNION_CFD_PROFILES.csv', 'auxiliary/robust_union/28_ROBUST_UNION_CFD_PROFILES.json', 'auxiliary/robust_union/29_ROBUST_UNION_EPSILON_COSTS.csv', 'auxiliary/robust_union/29_ROBUST_UNION_EPSILON_COSTS.json', 'auxiliary/robust_union/33_ROBUST_UNION_THRESHOLD_FRONTIERS.csv', 'auxiliary/robust_union/33_ROBUST_UNION_THRESHOLD_FRONTIERS.json', 'auxiliary/robust_union/35_ROBUST_UNION_EDGE_COST_SENSITIVITY.csv', 'auxiliary/robust_union/35_ROBUST_UNION_EDGE_COST_SENSITIVITY.json')

def selected_outputs(include_auxiliary=False):
    return MANUSCRIPT_DERIVED+(AUXILIARY_DERIVED if include_auxiliary else ())

def execute(directory,*args):
    subprocess.run([sys.executable,'-B',*args],cwd=directory,check=True)

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--work-dir',type=Path,
        help='New, nonexistent directory outside the extracted archive; retained for inspection.')
    parser.add_argument('--include-robust-union',action='store_true',help='Also reconstruct the separately labelled auxiliary robust-union analyses.')
    args=parser.parse_args()
    derived=selected_outputs(args.include_robust_union)
    execute(ROOT,'verify_manifest.py')
    if args.work_dir:
        work=args.work_dir.resolve()
        if work==ROOT or ROOT in work.parents:
            raise ValueError('The reconstruction directory must be outside the archive.')
        work.mkdir(parents=True,exist_ok=False)
    else:
        work=Path(tempfile.mkdtemp(prefix='pathway_reconstruction_')).resolve()
    print('Retained reconstruction workspace:',work,flush=True)
    # Copy only verified release contents. Never execute a writer in ROOT.
    with (ROOT/'SHA256SUMS.csv').open(newline='',encoding='utf-8-sig') as f:
        rows=list(csv.DictReader(f))
    for row in rows:
        target=work/row['path'];target.parent.mkdir(parents=True,exist_ok=True)
        shutil.copy2(ROOT/row['path'],target)
    shutil.copy2(ROOT/'SHA256SUMS.csv',work/'SHA256SUMS.csv')
    before={r['path']:r['sha256'] for r in rows}
    env=os.environ.copy();env['PYTHONDONTWRITEBYTECODE']='1'
    subprocess.run([sys.executable,'-B','verification_worker.py']+(['--include-robust-union'] if args.include_robust_union else []),cwd=work,env=env,check=True)
    comparisons={name:compare_file(ROOT/name,work/name) for name in derived}
    changed=[]
    for name,sha in before.items():
        current=hashlib.sha256((work/name).read_bytes()).hexdigest()
        if current!=sha:
            changed.append(name)
            if name not in derived:
                raise RuntimeError('Unexpected reconstructed-file change: '+name)
    report={'comparison_groups':{'manuscript':list(MANUSCRIPT_DERIVED),'auxiliary_robust_union':list(AUXILIARY_DERIVED) if args.include_robust_union else []},'original_archive_unchanged':None,'relative_tolerance':1e-10,
            'absolute_tolerance':1e-14,'comparisons':comparisons,
            'changed_bytes_in_work_copy':changed,
            'boundary_policy':'Discrete witnesses, counts, labels, classifications and certificate costs must agree. Boundary-sensitive disagreement is reported as failure, not relaxed.',
            'version':'2.1.9'}
    execute(ROOT,'verify_manifest.py')
    report['original_archive_unchanged']=True
    (work/'RECONSTRUCTION_REPORT.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Pathway-dependence certificates v2.1.9 verification passed',flush=True)

if __name__=='__main__':
    main()
