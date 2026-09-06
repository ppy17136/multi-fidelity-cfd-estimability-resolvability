"""Enumerate every connected partition of released CFD vectors before varying costs."""
from pathlib import Path
import hashlib,itertools,json,sys
import numpy as np
from weight_stability import envelope_bounds

HERE=Path(__file__).resolve().parent
DATA=HERE
sys.path.insert(0,str(DATA))
from run_real_cfd_common_partition_profiles import canonical_partition,allowed_indices

def main():
    output=[]; hashes={}
    for path in sorted((DATA/'data/cfd_mapping').glob('*.npz')):
        hashes[path.name]=hashlib.sha256(path.read_bytes()).hexdigest()
        with np.load(path) as x:
            raw=x['weighted_cell_vectors']; p=len(x['pathways']);q=len(x['fidelities'])
            names=[str(a) for a in x['fidelities']]; co=x['coefficients']
            case=str(x['case_id'][0]);support=str(x['support'][0])
        ix={n:i for i,n in enumerate(names)}
        edges=[(ix['F000'],ix['F100']),(ix['F010'],ix['F110']),(ix['F000'],ix['F010']),(ix['F100'],ix['F110'])]
        choices=list(itertools.product(range(p),repeat=q)); cmap={v:i for i,v in enumerate(choices)}
        v=np.stack([sum(co[i]*raw[c[i]*q+i] for i in range(q)) for c in choices])
        # Direct pair differences avoid Gram cancellation in near-attainment tests.
        dist=np.zeros((len(v),len(v)))
        for j in range(len(v)):
            dv=v[j:]-v[j];dist[j,j:]=np.einsum('ij,ij->i',dv,dv)
        dist=dist+dist.T;Dind=float(np.sqrt(dist.max()))
        parts={canonical_partition(q,[e for j,e in enumerate(edges) if not mask&(1<<j)]) for mask in range(16)}
        assert len(parts)==12
        records=[]
        for part in sorted(parts):
            ids=allowed_indices(part,p,cmap);diam=float(np.sqrt(dist[np.ix_(ids,ids)].max()))
            block={i:k for k,g in enumerate(part) for i in g}
            cut=[int(block[i]!=block[j]) for i,j in edges]
            records.append({'partition':part,'cut':cut,'diameter':diam,'relative_gap':max(0.,(Dind-diam)/Dind),'qualifies':diam>=(1-1e-7)*Dind})
        qualified=[r for r in records if r['qualifies']]
        # Remove componentwise dominated cuts only after checking all twelve partitions.
        pareto=[r for r in qualified if not any(all(a<=b for a,b in zip(t['cut'],r['cut'])) and t['cut']!=r['cut'] for t in qualified)]
        bounds=envelope_bounds([r['cut'] for r in pareto],[.15]*4,[.35]*4)
        output.append({'case_id':case,'support':support,'edge_order':['mesh_bottom','mesh_top','closure_left','closure_right'],
            'partitions':records,'pareto_cuts':[r['cut'] for r in pareto],'bounds':bounds})
        print(case,support,'cuts',len(pareto),bounds['minimum'],bounds['maximum_reoptimized'],bounds['minimum_fixed_witness_worst_cost'],flush=True)
    assert len(output)==8
    result={'scope':'fixed released CFD vectors; relative diameter tolerance 1e-7; no new CFD',
        'weight_domain':{'lower':[.15]*4,'upper':[.35]*4,'sum':1,'choice':'declared sensitivity domain, not estimated uncertainty'},
        'input_sha256':hashes,'audit_count':len(output),'results':output}
    (HERE/'weight_stability_results.json').write_bytes(json.dumps(result,indent=2).encode('utf-8'))

if __name__=='__main__':main()
