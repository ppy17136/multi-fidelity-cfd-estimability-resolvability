"""Evaluate all declared masks; do not select outputs by effect size."""
import os
for k in ('OMP_NUM_THREADS','OPENBLAS_NUM_THREADS','MKL_NUM_THREADS'):
    os.environ[k]='1'
from pathlib import Path
import csv, hashlib, json, sys, zipfile
import numpy as np
HERE=Path(__file__).resolve().parent
PACKAGE=HERE.parent
sys.path.insert(0,str(PACKAGE))
from graph_constrained_coupling import pathway_pair_compatibility, tree_graph_relaxation

def dump(p,x):p.write_bytes(json.dumps(x,indent=2,allow_nan=False).encode('utf-8'))
def writecsv(p,rows):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def groups(mask):
    starts=[0]+[i+1 for i in range(7) if mask&(1<<i)]
    return [list(range(a,b)) for a,b in zip(starts,starts[1:]+[8])]
def main():
    d=HERE/'results'
    gen=json.loads((d/'generation_summary.json').read_bytes())
    for n,key in [('raw_qoi.npz','raw_qoi_sha256'),('cellwise_contributions.npz','cellwise_sha256')]:
        assert hashlib.sha256((d/n).read_bytes()).hexdigest()==gen[key]
    raw=np.load(d/'raw_qoi.npz')['qoi']; arr=np.load(d/'cellwise_contributions.npz')
    Y=arr['Y'];u=arr['direction'];z=Y@u
    # Reconstruct via a different order: pathway means first, then differences.
    means=raw.mean(axis=1)
    check=np.zeros_like(Y)
    check[0]=means[:,0]
    for l in range(1,8):check[l]=means[:,l]-means[:,l-1]
    assert np.max(abs(check-Y))<1e-15
    allowed=pathway_pair_compatibility(z)
    edges=[(i,i+1,1/7) for i in range(7)]
    dp=tree_graph_relaxation(allowed,edges)
    cuts=[i for i in range(7) if dp.labels[i]!=dp.labels[i+1]]
    dpmask=sum(1<<i for i in cuts)
    masks=[];max_enum_error=0.
    for mask in range(128):
        bs=groups(mask)
        blocks=[z[b].sum(axis=0) for b in bs]
        lower=float(sum(v.min() for v in blocks));upper=float(sum(v.max() for v in blocks))
        values=np.array([0.])
        for v in blocks:values=(values[:,None]+v[None,:]).ravel()
        error=max(abs(lower-values.min()),abs(upper-values.max()))
        max_enum_error=max(max_enum_error,float(error))
        assert error<1e-13
        # Structural feasibility is based on exact stored-array extremizer
        # intersections, not an arbitrary tolerance on a summed width.
        supported=all(bool(set.intersection(*(set(allowed[i]) for i in b))) for b in bs)
        masks.append(dict(mask=mask,cuts=mask.bit_count(),cost=mask.bit_count()/7,
             lower=lower,upper=upper,width=upper-lower,full_width_supported=supported,
             components=';'.join(','.join(map(str,b)) for b in bs)))
    optimum=min(r['cost'] for r in masks if r['full_width_supported'])
    assert abs(dp.cost-optimum)<1e-14
    independent=float(np.ptp(z,axis=1).sum())
    coherent=float(np.ptp(z.sum(axis=0)))
    witness=masks[dpmask]
    assert abs(witness['width']-independent)<1e-13 and witness['full_width_supported']
    budget=[]
    for k in range(8):
        candidates=[r for r in masks if r['cuts']<=k]
        wr=max(candidates,key=lambda r:r['width'])
        budget.append(dict(cuts_budget=k,budget=k/7,
            maximum_common_partition_width=wr['width'],witness_mask=wr['mask'],
            lowest_endpoint=min(r['lower'] for r in candidates),
            highest_endpoint=max(r['upper'] for r in candidates)))
    def threshold(t):
        lo=masks[0]['lower'];hi=masks[0]['upper']
        state='above' if lo>t else 'below' if hi<t else 'indeterminate'
        if state=='indeterminate':return state,0.,0
        candidates=[r for r in masks if (r['lower']<=t if state=='above' else r['upper']>=t)]
        if not candidates:return state,None,None
        r=min(candidates,key=lambda r:r['cost'])
        return state,r['cost'],r['mask']
    points=sorted(set(r[k] for r in masks for k in ['lower','upper']))
    profiles=[]
    for j,t in enumerate(points):
        state,cost,mask=threshold(t)
        profiles.append(dict(kind='point',left=t,right=t,baseline_state=state,
                             minimum_cost=cost,witness_mask=mask))
        if j<len(points)-1:
            right=points[j+1];mid=t+(right-t)/2
            # Do not invent a representable interior point in an adjacent-float gap.
            if not t<mid<right:continue
            state,cost,mask=threshold(mid)
            profiles.append(dict(kind='open_interval',left=t,right=right,baseline_state=state,
                                 minimum_cost=cost,witness_mask=mask))
    summary=dict(cells=8,pathways=4,output_dimension=4,graph='seven-edge chain',
        masks_checked=128,independent_assignments=65536,
        direction=u.tolist(),coherent_min=masks[0]['lower'],coherent_max=masks[0]['upper'],
        independent_min=masks[-1]['lower'],independent_max=masks[-1]['upper'],
        coherent_width=coherent,independent_width=independent,
        width_expansion=independent/coherent if coherent else None,
        directional_cost=dp.cost,exhaustive_mask_cost=optimum,
        cut_edge_indices=cuts,labels=[list(x) for x in dp.labels],
        label_list_sizes=[len(x) for x in allowed],
        label_lists=[[list(x) for x in a] for a in allowed],
        all_labels_singleton=all(len(x)==1 for x in allowed),
        witness_width=witness['width'],endpoint_enumeration_max_error=max_enum_error,
        decomposition_reconstruction_max_error=float(np.max(abs(check-Y))),
        above_baseline_fragility_interval=[masks[-1]['lower'],masks[0]['lower']],
        below_baseline_fragility_interval=[masks[0]['upper'],masks[-1]['upper']],
        outside_independent_interval='no crossing possible; minimum cost infinite',
        already_indeterminate='minimum cost zero; not a newly destroyed decision',
        note='All frozen outcomes retained. Singleton lists imply a predetermined cut, not a difficult DP optimization.',
        graph_module_sha256=hashlib.sha256((PACKAGE/'graph_constrained_coupling.py').read_bytes()).hexdigest(),
        vector_diameter_computed=False)
    expected=json.loads((d/'certificate_summary.json').read_bytes())
    for key,value in summary.items():
        if isinstance(value,float):
            assert abs(value-expected[key])<1e-13,(key,value,expected[key])
        else:
            assert value==expected[key],key
    for filename,actual in [('all_partition_profiles.csv',masks),
                            ('budget_profiles.csv',budget),
                            ('threshold_profiles.csv',profiles)]:
        with (d/filename).open(newline='',encoding='utf-8') as f:
            recorded=list(csv.DictReader(f))
        assert len(actual)==len(recorded),filename
        for a,b in zip(actual,recorded):
            for key,value in a.items():
                if value is None:
                    assert b[key]==''
                elif isinstance(value,float):
                    assert abs(value-float(b[key]))<1e-13,(filename,key)
                else:
                    assert str(value)==b[key],(filename,key)
    freeze=json.loads((d/'FREEZE.json').read_bytes())
    with zipfile.ZipFile(PACKAGE/'provenance/diffusion_precomputation.zip') as original:
        assert set(original.namelist())=={'PROTOCOL.md','run_benchmark.py'}
        for name,key in [('PROTOCOL.md','protocol_sha256'),('run_benchmark.py','generator_sha256')]:
            assert hashlib.sha256(original.read(name)).hexdigest()==freeze[key]
    assert hashlib.sha256((d/'pathway_nodes.npz').read_bytes()).hexdigest()==freeze['nodes_sha256']
    print('Stochastic diffusion: frozen inputs, 128 partitions, budget and threshold profiles verified.')
if __name__=='__main__':main()
