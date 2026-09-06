"""Check the CSV-backed manuscript verification table and test inventory.

This check does not rerun the frozen timing benchmarks. run_all.py executes
the unit tests before invoking it; the table check separately verifies their
discovered counts and the rendered CSV rows.
"""
from pathlib import Path
import argparse
import csv
import re
import unittest

ROOT=Path(__file__).resolve().parent

def read_rows():
    with (ROOT/'45_ALGORITHM_CERTIFICATE_SUMMARY.csv').open(newline='',encoding='utf-8') as f:
        return list(csv.DictReader(f))

def render(rows):
    header = r"""\begin{tabularx}{\textwidth}{@{}p{0.19\textwidth}r p{0.24\textwidth}X@{}}
\toprule
Verification block & Instances & Frozen scale & Result \\
\midrule
"""
    body = ''.join(' & '.join([r['verification_block'],format(int(r['instances']),','),r['display_scale_tex'],r['display_result_tex']])+r' \\'+'\n' for r in rows)
    return header + body + '\\bottomrule\n\\end{tabularx}\n'

def validate_rows(rows,total,weight):
    names=[r['verification_block'] for r in rows]
    if len(rows)!=8 or len(set(names))!=8:
        raise ValueError('Expected eight distinct verification blocks.')
    by_name={r['verification_block']:r for r in rows}
    for name,count in [('Core regression suite',total-weight),('Weight stability',weight)]:
        r=by_name[name]
        if int(r['instances'])!=count:
            raise ValueError(f'{name}: CSV count differs from discovered tests.')
        for key in ['result','display_result_tex']:
            m=re.match(r'(\d+)/(\d+)',r[key].replace(',', ''))
            if not m or tuple(map(int,m.groups()))!=(count,count):
                raise ValueError(f'{name}: inconsistent passing-count statement.')
    for r in rows:
        for key in ['result','display_result_tex']:
            m=re.match(r'(\d+)/(\d+)',r[key].replace(',', ''))
            if not m or int(m[2])!=int(r['instances']) or int(m[1])!=int(m[2]):
                raise ValueError(f"{r['verification_block']}: inconsistent summary count.")

def main():
    ap=argparse.ArgumentParser();ap.add_argument('--write',action='store_true');args=ap.parse_args()
    loader=unittest.TestLoader()
    total=loader.discover(str(ROOT),pattern='test_*.py').countTestCases()
    weight=loader.discover(str(ROOT),pattern='test_weight_stability.py').countTestCases()
    if loader.errors:raise RuntimeError('\n'.join(loader.errors))
    rows=read_rows();validate_rows(rows,total,weight)
    text=render(rows);p=ROOT/'45_ALGORITHM_CERTIFICATE_TABLE_ROWS.tex'
    if args.write:p.write_text(text,encoding='utf-8')
    if not p.is_file() or p.read_text(encoding='utf-8')!=text:
        raise RuntimeError('CSV and manuscript table fragment differ.')
    print(f'Algorithm table consistency passed: {len(rows)} rows; {total-weight} core + {weight} weight tests.')

if __name__=='__main__':main()
