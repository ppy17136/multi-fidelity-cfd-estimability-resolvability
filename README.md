# From predictive accuracy to decision-ready multi-fidelity simulation

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21670992.svg)](https://doi.org/10.5281/zenodo.21670992)

This repository contains the reproducibility materials for:

> *From predictive accuracy to decision-ready multi-fidelity simulation: Minimum-cost contrast-support repair and dependence-aware decisions*

## Central contribution

The method treats scientific attribution as a sequential support-and-resolution problem rather than as predictive validation alone:

1. verify that each declared contrast is estimable from the acquired design;
2. repair missing support using exact finite-pool minimum-cost enumeration or an explicitly approximate greedy solver;
3. propagate contrast-aligned numerical uncertainty under a declared dependence construction; and
4. acquire the next computation until the claim is resolved above the minimum, resolved below the minimum, or retained as indeterminate.

Release `v1.3.0` corrects the closed-loop comparison by removing a duplicate policy alias and assigning all seven distinct policies the same action-indexed random scenario within each truth-replicate pair. It also retains the exhaustive complete-lattice dependence audit: four inventory-marked complete periodic-hill geometries on two supports, giving eight audits. Because the coherent set is nested in the independently recombined set, non-increasing robustness ratio is structural; the empirical findings are the 10.83--15.71-fold exact finite-set diameter expansion and seven paired mapping-robustness classification reversals at the frozen `R = 3` rule.

## Latest reproducibility package

The latest self-contained package is in `v20_method_upgrade/`.

```bash
python -m pip install -r v20_method_upgrade/requirements.txt
cd v20_method_upgrade
python run_all.py
```

A clean Linux extraction runs 17 automated tests, reconstructs the original four-audit result, and independently rebuilds all eight complete-lattice audits from released float64 weighted-vector caches. The maximum absolute difference from expected output is approximately `2.66e-14`.

## Reproducibility boundary

The released vectors reproduce the finite-set mapping-robustness classifications exactly. They do not independently rebuild omitted raw OpenFOAM cache eligibility, rerun the underlying OpenFOAM simulations, or represent total CFD uncertainty. The exact minimum-cost certificate applies only to Gate-1 repair over the declared finite candidate pool; greedy repair and total cost to a terminal decision are not globally certified.

## Repository map

- `v20_method_upgrade/`: current V20 protocols, code, corrected paired-policy outputs, tests, eight float64 audit caches, expected outputs, figures, and hashes.
- `v16_method_upgrade/`: preceding complete-lattice release retained for provenance.
- `v9_method_upgrade/`: earlier method-upgrade release retained for provenance.
- `theory/`, `protocols/`, `code/`, `data/`, `figures/`, and `provenance/`: earlier ordered two-gate materials.
- `THIRD_PARTY_DATA_NOTICE.md`: data ownership and redistribution boundaries.

## Data boundaries

This archive contains authored code, protocols, compact derived CFD audit data, figures, hashes, and completion evidence. It does not contain downloaded articles, third-party DNS payloads, raw processor directories, or full CFD time histories.

## Citation and DOI

Repository: https://github.com/ppy17136/multi-fidelity-cfd-estimability-resolvability

- Concept DOI for all versions: https://doi.org/10.5281/zenodo.21670992
- Each GitHub release is archived by Zenodo under a distinct immutable version DOI.

## Licences

- Source code: MIT (`LICENSE-CODE`).
- Authored data, figures, and documentation: CC BY 4.0 (`LICENSE-DATA`).
- Third-party data are not redistributed.
