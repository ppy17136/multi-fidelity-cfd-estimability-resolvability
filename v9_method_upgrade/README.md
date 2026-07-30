# V9 claim-adaptive method upgrade

This directory contains the frozen method-upgrade protocol, exact and greedy
support-repair implementations, deterministic and probabilistic resolution
rules, closed-loop acquisition benchmarks, real-CFD mapping audits, a
second-physics transfer audit, machine-readable outputs, and tests.

Key verified outcomes:

- exact support repair restored every declared contrast in six benchmarks;
- greedy repair was feasible in all six and matched exact cost in five;
- bounded-error Gate-2 calibration produced no false decisions in 200,000
  realisations per scenario;
- the proposed state switch reduced median valid-decision cost by 35% in the
  frozen subthreshold 3x2 benchmark;
- four real-CFD mapping cases passed under a coherent joint ensemble and failed
  under an independent-cell envelope; and
- the pressure-loss transfer is estimable but remains explicitly unresolved
  because no independent contrast-level error bound was available.

Run `python -m pytest -q` for the 12-test verification suite. See
`REPRODUCIBILITY_BOUNDARY.md` for raw-data and scope limitations.
