# Frozen second-physics transfer audit

Protocol: `V9-20260730`

## Purpose

Test whether the support/contrast part of the proposed method transfers from
periodic-hill velocity fields to an already completed obstruction-flow CFD
family with a scalar response and different fidelity axes. This is a
retrospective transfer audit, not a new prospective CFD experiment.

## Frozen design

- Physical class: circular internal flow with inlet or outlet obstruction.
- Conditions:
  - unobstructed reference (`normal`);
  - severe inlet obstruction (`severe_inlet`);
  - severe outlet obstruction (`severe_outlet`).
- Fidelity axis 1: mesh, `medium` versus `fine`.
- Fidelity axis 2: turbulence closure, `kEpsilon` versus `kOmegaSST`.
- Response: boundary pressure-loss coefficient `K_boundary`.
- Contrast:

  \[
  C_K =
  (K_{\mathrm{fine,SST}}-K_{\mathrm{fine},k\varepsilon})
  -
  (K_{\mathrm{medium,SST}}-K_{\mathrm{medium},k\varepsilon}).
  \]

## Eligibility

The audit is included only if all four cell means are present for every
condition and each archived calculation is marked as ended with no material
mass imbalance.

## Claim boundary

- The audit can establish design-supported estimability and quantify the
  mesh-by-closure contrast.
- It cannot establish a risk-controlled Gate-2 decision because the archive
  does not provide an independent, contrast-level bound on pressure-loss
  error for all four cells.
- Solver residuals and mass imbalance are numerical-quality diagnostics; they
  are not converted into uncertainty in `K_boundary`.
- The prior manuscript containing related pump data was rejected and is not a
  published source. Any later reuse must still be disclosed internally and
  checked for repository/version overlap.
