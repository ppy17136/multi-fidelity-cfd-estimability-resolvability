# Frozen real-CFD contrast-mapping repair implementation

Date frozen: 2026-07-30

This implementation applies to the two Protocol-179 primary geometries:

- `alph15-7929-2024`;
- `alph15-13929-2024`.

The archived periodic-Delaunay linear map remains pathway A. Pathway B is an
independently implemented periodic inverse-distance-weighted nearest-neighbour
map in the archived transformed coordinates \((x/L_x,\eta)\).

Primary pathway-B settings:

- periodic copies at \(x/L_x-1\), \(x/L_x\), and \(x/L_x+1\);
- \(k=8\) neighbours;
- inverse-distance exponent \(p=2\);
- exact-coordinate matches use the matched value without distance weighting;
- the common query support and physical-volume weights are the finite
  pathway-A F100 support used by the archived four-cell contrast.

For each fidelity cell \(i\), the pathway sensitivity is

\[
e_i = \left\|u_i^{(A)}-u_i^{(B)}\right\|_{L^2_V}.
\]

For the mixed contrast with coefficients \((+1,-1,-1,+1)\), the propagated
pathway-sensitivity bound is

\[
B_C=\sum_i |c_i|e_i.
\]

This is a deterministic bound on the difference between the two implemented
mapping pathways by the triangle inequality. It is not a bound on total CFD
error, model discrepancy, or distance to the unknown physical truth.

Robustness-only settings are \(k=4\) and \(k=16\), with all other choices
unchanged. They do not replace the frozen primary \(k=8\) result.
