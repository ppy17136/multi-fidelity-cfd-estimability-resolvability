"""Independently enumerate the three-cell example without plotting libraries."""
from pathlib import Path
import itertools
import json

ROOT=Path(__file__).resolve().parent

def main():
    frozen=json.loads((ROOT/'worked_example_verified.json').read_text())
    z=frozen['contributions']
    partitions=[[[0,1,2]],[[0],[1,2]],[[0,1],[2]],[[0],[1],[2]]]
    records=[]
    costs=[0,.5,.5,1]
    for p in partitions:
        values=sorted({sum(z[i][k] for group,k in zip(p,ks) for i in group)
                       for ks in itertools.product(range(3),repeat=len(p))})
        records.append({'partition':p,'values':values,'lower':min(values),
                        'upper':max(values),'width':max(values)-min(values)})
    assert records==frozen['partitions']
    assert [r['width'] for r in records]==[10,11,13,14]
    assert max(records[i]['width'] for i in [0,1,2])==13
    assert max(records[i]['upper'] for i in [0,1,2])-min(records[i]['lower'] for i in [0,1,2])==14
    # The scalar interval rule retains equality as indeterminate.
    eta=min(c for c,r in zip(costs,records) if r['lower']<=-2<=r['upper'])
    rho=min(c for c,r in zip(costs,records) if r['width']==records[-1]['width'])
    assert eta==frozen['eta']==.5 and rho==frozen['rho']==1
    print('Worked-example verification passed: four partitions; common width 13 versus union width 14; eta=0.5, rho=1.')

if __name__=='__main__':main()
