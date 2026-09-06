# Stochastic diffusion benchmark specification

This is a new computational benchmark, not an experiment or an observed physical dataset. The reported target is directional width, not full vector diameter.

This explanatory edition describes the released numerical specification.
The original precomputation protocol and generator are retained byte-for-byte
in provenance/diffusion_precomputation.zip and checked against the original
FREEZE.json. This edition does not replace that frozen record. The following
specification was recorded before generation of the response array. All results, including zero-width gaps, zero/full-cost certificates and no decision changes, will be retained. No selection by effect size is allowed.

## Model

Domain [0,1]^2, homogeneous Dirichlet boundary, forcing 1:

    -div(a(x,xi) grad(v)) = 1.
    a = exp(0.7 * sum_{k=1}^5 k^(-3/2) xi_k phi_k(x)).
    xi_k independent Uniform[-1,1].

Modes, in order:

    sin(pi*x) sin(pi*y)
    cos(2*pi*x) sin(pi*y)
    sin(pi*x) cos(2*pi*y)
    sin(2*pi*x) sin(2*pi*y)
    cos(3*pi*x) cos(pi*y)

Eight uniform Cartesian finite-volume grids have N=8,12,16,24,32,40,48,64 cells per side. Cell-centred unknowns; harmonic interior-face conductivity; positive analytic boundary-face conductivity and half-cell Dirichlet distance; double precision sparse direct solve. No estimated runtime is treated as a scientific observation.

## Frozen numerical pathways

Four equally weighted 32-point integration rules in five dimensions:

- scrambled Sobol, seed 2026090601, base-two generation m=5;
- scrambled Halton, seed 2026090602;
- Latin hypercube, seed 2026090603;
- tensor two-point Gauss-Legendre, nodes +/-1/sqrt(3).

Nodes are shared across all mesh levels within a pathway and saved before PDE solves. Each rule has weight 1/32. These are frozen alternative integration rules, not statistical replicates or a confidence distribution.

## Outputs and decomposition

Four cell-volume-overlap weighted area means over:

    [0.20,0.35] x [0.20,0.35]
    [0.65,0.80] x [0.20,0.35]
    [0.20,0.35] x [0.65,0.80]
    [0.65,0.80] x [0.65,0.80].

Fractional overlaps keep the physical areas fixed across grids. This is a piecewise-constant reconstruction of each finite-volume solution, not an exact continuum integral.

    Y[0,p] = A_p Q_0;
    Y[l,p] = A_p(Q_l - Q_{l-1}), l>0.

Check sum_l Y[l,p] = A_p Q_finest. The graph is a chain of the eight signed correction cells; seven edges each cost 1/7. Retained edges enforce the declared common integration-rule identity on adjacent correction cells. This is not a universal requirement for valid MLMC estimators.

Primary direction u=(1,1,1,1)/2 is fixed. Use exact stored-array extrema (no tie enlargement). Report every extremizer-list size, returned witness, direct witness width, complete 128-mask lower/upper/width profiles, and full threshold intervals and equality boundaries. If each list is a singleton, report that the optimization is structurally trivial.

## Numerical checks and sequencing

1. Manufacture v=sin(pi*x)sin(pi*y), with a=1 and a=exp(0.3*x). Analytic forcing for the variable case is a*(2*pi^2*v - 0.3*pi*cos(pi*x)*sin(pi*y)). Check positive solution, relative residual <=1e-10, decreasing L2 errors over N=8,16,32,64 and last refinement order >1.5. Stop for investigation on failure; do not silently change the model.
2. Time eight frozen stochastic nodes (first two from each pathway) at every grid: 64 solves. Extrapolation from this pilot is labelled an estimate. The pilot uses a subset of the final nodes and is not an independent confirmation.
3. If correct and modest in cost, compute all 1024 solves locally, using one sequential CPU worker and single-thread native libraries. Save all QoIs, residuals and elapsed solve times.
4. Check telescoping relative error <=1e-12 and all finite values. Compare the existing forest implementation with exhaustive 128-mask checks based on exact label intersections. Check each partition's endpoints independently with scalar blockwise pathway enumeration (at most 4^8 values). Do not form an all-pairs vector-distance matrix.

The released responses include the completed numerical checks and independent
solver audit. Run python stochastic_diffusion/run_benchmark.py --work-dir
NEW_EXTERNAL_DIRECTORY --full from the package root for a new full generation.
The launcher copies the unchanged original sources to the new directory and
executes pilot checks followed by full generation. See the package README
for read-only verification of the released arrays.
