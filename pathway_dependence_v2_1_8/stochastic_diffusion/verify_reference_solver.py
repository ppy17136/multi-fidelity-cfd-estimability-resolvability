"""Read-only repetition of the independent implementation checks."""
from independent_solver_audit import reference, D
import numpy as np
nodes=np.load(D/'pathway_nodes.npz')['nodes']
raw=np.load(D/'raw_qoi.npz')['qoi']
maximum=0.
for p in range(4):
    for sample in (0,31):
        for level,n in ((0,8),(2,16),(7,64)):
            q,residual=reference(n,nodes[p,sample])
            difference=float(np.max(np.abs(q-raw[p,sample,level])))
            assert difference<1e-10
            maximum=max(maximum,difference)
print('Independent solver: 24 checks passed; max absolute QoI difference',maximum)
