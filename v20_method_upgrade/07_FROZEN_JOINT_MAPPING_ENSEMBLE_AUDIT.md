# Frozen joint mapping-ensemble uncertainty audit

Date frozen: 2026-07-30

The declared mapping ensemble contains exactly four prespecified pathways:

1. periodic physical-coordinate Delaunay-linear interpolation;
2. transformed-coordinate periodic IDW with \(k=4,p=2\);
3. transformed-coordinate periodic IDW with \(k=8,p=2\);
4. transformed-coordinate periodic IDW with \(k=16,p=2\).

The ensemble is evaluated on both previously declared supports:

- F100 finite common support;
- F000 DNS-eligible shared representable support.

For each pathway \(m\), the four-cell mixed contrast field is \(I_m\). The
finite-ensemble diameter is

\[
D_{\mathcal M}=\max_{m,n}\|I_m-I_n\|_{L^2_V}.
\]

The conservative ensemble signal is

\[
S_{\mathcal M}=\min_m\|I_m\|_{L^2_V}.
\]

The prespecified descriptive resolution ratio is
\(S_{\mathcal M}/D_{\mathcal M}\), with the inherited threshold 3.

This is a deterministic robustness statement only over the declared finite
mapping ensemble. It is not a universal interpolation-error bound, a
probabilistic confidence statement, or a bound on total CFD/model error. The
axis-aligned cell-wise triangle bound remains reported as a separate, more
adversarial uncertainty set.
