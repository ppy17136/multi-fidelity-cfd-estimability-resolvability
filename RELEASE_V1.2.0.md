# Release v1.2.0 - Complete-lattice dependence audit

This release adds the self-contained V16 complete-lattice reconstruction supporting the revised manuscript.

## New evidence

- Objective inclusion of every eligible archived four-cell lattice: four geometries x two supports = eight audits.
- At frozen `R = 3`: coherent sets pass in 7/8 audits; independent recombinations pass in 0/8; seven paired reversals.
- `R_recombined < R_coherent` in 8/8 audits.
- Exact finite-set diameter expansion: 10.8281-15.7142 times.
- Common diagnostic separation interval: `1.6442596248 < tau < 2.2831676147`.
- Eight released float64 weighted-vector caches and a frozen eligibility index.
- Clean-extraction one-command verification: 15 tests passed; maximum expected-output difference approximately `2.66e-14`.

## Reproduce

```bash
cd v16_method_upgrade
python run_all.py
```

## Boundary

The package reconstructs the numerical decision from released weighted vectors. It does not rerun the underlying OpenFOAM simulations and does not claim to represent total CFD uncertainty.

Concept DOI: https://doi.org/10.5281/zenodo.21670992