"""Independent scalar-loop assembly and iterative-solver audit of frozen data.

Shares the declared discretization, not the generator implementation. This is
implementation verification, not independent physical validation.
"""
import os
for name in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS'):
    os.environ[name] = '1'
from pathlib import Path
import hashlib
import json
import math
import time
import numpy as np
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import cg

ROOT = Path(__file__).resolve().parent
D = ROOT / 'results'

def conductivity(x, y, xi):
    p = math.pi
    phi = (math.sin(p*x)*math.sin(p*y),
           math.cos(2*p*x)*math.sin(p*y),
           math.sin(p*x)*math.cos(2*p*y),
           math.sin(2*p*x)*math.sin(2*p*y),
           math.cos(3*p*x)*math.cos(p*y))
    return math.exp(0.7*sum(float(xi[k])*phi[k]/(k+1)**1.5 for k in range(5)))

def reference(n, xi):
    h = 1.0/n
    a = np.array([[conductivity((i+.5)*h,(j+.5)*h,xi)
                   for i in range(n)] for j in range(n)])
    matrix = lil_matrix((n*n, n*n))
    for j in range(n):
        for i in range(n):
            row = j*n+i
            for di,dj in ((-1,0),(1,0),(0,-1),(0,1)):
                ii,jj = i+di,j+dj
                if 0 <= ii < n and 0 <= jj < n:
                    conductance = 2*a[j,i]*a[jj,ii]/(a[j,i]+a[jj,ii])/h**2
                    matrix[row,row] += conductance
                    matrix[row,jj*n+ii] -= conductance
                else:
                    x = (i+.5+.5*di)*h
                    y = (j+.5+.5*dj)*h
                    matrix[row,row] += 2*conductivity(x,y,xi)/h**2
    matrix = matrix.tocsr()
    rhs = np.ones(n*n)
    solution, info = cg(matrix, rhs, rtol=1e-12, atol=0, maxiter=20*n*n)
    residual = float(np.linalg.norm(matrix@solution-rhs)/np.linalg.norm(rhs))
    assert info == 0 and residual < 2e-11
    q = []
    for x0,x1,y0,y1 in ((.2,.35,.2,.35),(.65,.8,.2,.35),
                        (.2,.35,.65,.8),(.65,.8,.65,.8)):
        total = 0.0
        for j in range(n):
            for i in range(n):
                dx = max(0,min((i+1)*h,x1)-max(i*h,x0))
                dy = max(0,min((j+1)*h,y1)-max(j*h,y0))
                total += dx*dy*solution[j*n+i]
        q.append(total/((x1-x0)*(y1-y0)))
    return np.array(q), residual

def main():
    start = time.perf_counter()
    nodes = np.load(D/'pathway_nodes.npz')['nodes']
    raw = np.load(D/'raw_qoi.npz')['qoi']
    rows = []
    # Fixed first/last nodes of every pathway and coarse/intermediate/fine mesh.
    for p in range(4):
        for sample in (0,31):
            for level,n in ((0,8),(2,16),(7,64)):
                q,residual = reference(n,nodes[p,sample])
                difference = float(np.max(np.abs(q-raw[p,sample,level])))
                assert difference < 1e-10
                rows.append(dict(pathway=p,sample=sample,level=level,n=n,
                                 max_absolute_qoi_difference=difference,
                                 relative_residual=residual))
    report = dict(passed=True,checks=len(rows),rows=rows,
                  maximum_absolute_qoi_difference=max(r['max_absolute_qoi_difference'] for r in rows),
                  maximum_relative_residual=max(r['relative_residual'] for r in rows),
                  elapsed_s=time.perf_counter()-start,
                  script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  scope='Independent implementation and linear solver; same mathematical discretization; not physical validation')
    target = D/'independent_solver_audit.json'
    assert not target.exists(), 'Preserve existing audit; do not overwrite'
    target.write_text(json.dumps(report,indent=2,allow_nan=False),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k!='rows'},indent=2))

if __name__ == '__main__':
    main()
