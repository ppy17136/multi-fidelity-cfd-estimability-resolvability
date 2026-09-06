"""Numerical agreement is separate from the frozen archive's byte identity."""
import csv
import json
import math
import re
from pathlib import Path

RTOL=1e-10
ATOL=1e-14
NUMBER=re.compile(r'^[+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?$')
DISCRETE=('cost','budget','epsilon','count','members','partitions','partition',
          'cut','witness','qualif','index','label','boundar','threshold_state')

def compare_values(expected,actual,path='',stats=None,strict=False):
    if stats is None:
        stats={'numeric_values':0,'maximum_absolute_difference':0.0}
    if isinstance(expected,dict):
        if not isinstance(actual,dict) or expected.keys()!=actual.keys():
            raise AssertionError('Object keys differ: '+path)
        for key in expected:
            compare_values(expected[key],actual[key],path+'/'+key,stats,
                           strict or key in ('partition','cut','pareto_cuts','minimum_diameter_cost_by_epsilon') or (not isinstance(expected[key],(dict,list)) and any(t in key.lower() for t in DISCRETE)))
    elif isinstance(expected,list):
        if not isinstance(actual,list) or len(expected)!=len(actual):
            raise AssertionError('List shape differs: '+path)
        for i,(a,b) in enumerate(zip(expected,actual)):
            compare_values(a,b,path+'/'+str(i),stats,strict)
    elif isinstance(expected,bool) or expected is None:
        if type(expected)!=type(actual) or expected!=actual:
            raise AssertionError('Logical value differs: '+path)
    elif isinstance(expected,(int,float)) and isinstance(actual,(int,float)) and not isinstance(actual,bool):
        stats['numeric_values']+=1
        if not math.isfinite(expected) or not math.isfinite(actual):
            if expected!=actual:raise AssertionError('Nonfinite value differs: '+path)
        elif strict or isinstance(expected,int):
            if expected!=actual:raise AssertionError('Discrete/certificate value differs: '+path)
        else:
            error=abs(expected-actual)
            stats['maximum_absolute_difference']=max(stats['maximum_absolute_difference'],error)
            if not math.isclose(expected,actual,rel_tol=RTOL,abs_tol=ATOL):
                raise AssertionError(f'Numeric mismatch {path}: {expected} vs {actual}')
    elif expected!=actual or type(expected)!=type(actual):
        raise AssertionError('Text/type differs: '+path)
    return stats

def csv_scalar(value):
    if NUMBER.fullmatch(value):
        return float(value) if '.' in value or 'e' in value.lower() else int(value)
    return value

def compare_file(frozen,recomputed):
    frozen=Path(frozen);recomputed=Path(recomputed)
    if frozen.suffix=='.json':
        a=json.loads(frozen.read_text(encoding='utf-8-sig'))
        b=json.loads(recomputed.read_text(encoding='utf-8-sig'))
    elif frozen.suffix=='.csv':
        def read(path):
            with path.open(newline='',encoding='utf-8-sig') as f:
                return [{k:csv_scalar(v) for k,v in row.items()} for row in csv.DictReader(f)]
        a,b=read(frozen),read(recomputed)
    else:
        raise ValueError('Only JSON/CSV numeric reconstruction is supported')
    return compare_values(a,b,frozen.name)
