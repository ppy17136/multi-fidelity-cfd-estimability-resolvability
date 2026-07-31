# Frozen shared-representable-support audit

Date frozen: 2026-07-30

The audit tests whether the four-cell fidelity-source contrast changes its
numerical interpretation when it is evaluated on the coarsest support shared
by all four cells, instead of on the F100 medium-grid support.

Frozen choices:

- geometries: `alph15-7929-2024` and `alph15-13929-2024`;
- common query support: DNS-eligible F000 physical cell centres;
- physical-volume weights: F000 cell volumes;
- pathway A: archived periodic Delaunay-linear interpolation;
- pathway B: periodic transformed-coordinate IDW with \(k=8,p=2\);
- robustness-only IDW settings: \(k=4\) and \(k=16\);
- mixed coefficients: \((+1,-1,-1,+1)\) for
  F000/F100/F010/F110.

The primary propagated pathway bound is the sum of the four cell-wise
pathway-difference RMS values. The direct difference between the two mixed
contrast fields is descriptive repeatability evidence only.

This audit does not redefine the original fine-support estimand after seeing
the result. Fine-support and shared-support contrasts will be reported
together. A favourable shared-support result permits attribution only in the
shared representable subspace; it does not validate resolution-exclusive
fine-scale content.
