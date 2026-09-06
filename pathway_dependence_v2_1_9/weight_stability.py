"""Finite-witness cost sensitivity; optimization is tolerance-qualified."""
import numpy as np
from scipy.optimize import linprog

def validate(cuts, lower, upper):
    d=np.asarray(cuts,dtype=float); lo=np.asarray(lower,dtype=float); hi=np.asarray(upper,dtype=float)
    if d.ndim!=2 or len(d)==0 or not np.isin(d,[0,1]).all():raise ValueError('nonempty binary witness matrix required')
    if lo.shape!=(d.shape[1],) or hi.shape!=lo.shape or not np.isfinite(lo).all() or not np.isfinite(hi).all():raise ValueError('invalid bounds')
    if np.any(lo<=0) or np.any(lo>hi) or lo.sum()>1+1e-12 or hi.sum()<1-1e-12:raise ValueError('infeasible strictly positive normalized weight bounds')
    return d,lo,hi

def lp(c,lo,hi,A=None,b=None):
    r=linprog(c,A_ub=A,b_ub=b,A_eq=np.ones((1,len(c))),b_eq=[1.],bounds=list(zip(lo,hi)),method='highs')
    if not r.success:raise RuntimeError(r.message)
    if abs(r.x.sum()-1)>1e-8 or np.min(r.x-lo)<-1e-8 or np.max(r.x-hi)>1e-8:raise RuntimeError('invalid LP witness')
    return float(r.fun),r.x

def envelope_bounds(cuts,lower,upper):
    d,lo,hi=validate(cuts,lower,upper)
    lows=[lp(a,lo,hi) for a in d]
    highs=[lp(-a,lo,hi) for a in d]
    # max_w min_a w.a is a hypograph LP, not min_a max_w w.a.
    A=np.column_stack([-d,np.ones(len(d))])
    r=linprog(np.r_[np.zeros(d.shape[1]),-1.],A_ub=A,b_ub=np.zeros(len(d)),
        A_eq=np.array([np.r_[np.ones(d.shape[1]),0.]]),b_eq=[1.],bounds=list(zip(lo,hi))+[(None,None)],method='highs')
    if not r.success:raise RuntimeError(r.message)
    reoptimized=float(np.min(d@r.x[:-1]))
    if abs(reoptimized-r.x[-1])>1e-8:raise RuntimeError('hypograph value mismatch')
    regions=[]
    for a in d:
        Areg=d*0+a-d
        feasible=linprog(np.zeros(d.shape[1]),A_ub=Areg,b_ub=np.zeros(len(d)),A_eq=np.ones((1,d.shape[1])),b_eq=[1.],bounds=list(zip(lo,hi)),method='highs')
        regret=max(-lp(-(a-b),lo,hi)[0] for b in d)
        regions.append({'cut':a.astype(int).tolist(),'optimal_somewhere':bool(feasible.success),'worst_regret':max(0.,regret)})
    fixed=min(-r[0] for r in highs)
    return {'minimum':min(r[0] for r in lows),'maximum_reoptimized':reoptimized,
        'minimum_fixed_witness_worst_cost':fixed,'commitment_gap':fixed-reoptimized,
        'maximizing_weights':r.x[:-1].tolist(),'witness_regions':regions}
