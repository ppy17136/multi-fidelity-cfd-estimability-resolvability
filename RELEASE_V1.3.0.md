# Release v1.3.0 - Corrected paired-policy benchmark

This release adds the self-contained V20 verification package supporting the associated study.

## Corrections and new evidence

- Removed a duplicate closed-loop policy alias; seven distinct policies are now reported.
- Assigned every policy the same pre-generated action-indexed random scenario within each truth-replicate pair.
- Regenerated the 5,000-replicate two-truth benchmark and the 2,000-replicate, 11-point truth-grid audit.
- Added automated checks for policy uniqueness and shared-scenario determinism.
- Retained the complete-lattice audit: eight audits, 7/8 coherent passes, 0/8 independently recombined passes, and seven paired mapping-robustness classification reversals at the frozen `R = 3` rule.
- Retained the 10.83--15.71-fold exact finite-set diameter expansion and the common diagnostic separation interval `1.6442596248 < tau < 2.2831676147`.
- Synchronized Figures 1--3 and Supplementary Figure S1 with the archived results.
- Normalized all ZIP member paths for cross-platform extraction.
- Clean Linux extraction: 17 tests passed; maximum expected-output difference approximately `2.66e-14`.

## Reproduce

```bash
cd v20_method_upgrade
python run_all.py
```

## Boundaries

The exact minimum-cost certificate applies only to Gate-1 support restoration over the declared finite pool. Greedy repair is approximate, and the closed-loop policy is not claimed to minimize total cost globally. The released weighted vectors reproduce finite-set mapping-robustness classifications; they do not independently rebuild raw OpenFOAM cache eligibility or total CFD uncertainty.

Concept DOI: https://doi.org/10.5281/zenodo.21670992
