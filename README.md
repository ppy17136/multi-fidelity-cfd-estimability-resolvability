# Contrast-directed multi-fidelity design for CFD

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21670992.svg)](https://doi.org/10.5281/zenodo.21670992)

This repository contains the reproducibility materials for:

> *Contrast-directed multi-fidelity design for CFD: Minimum-cost support repair and dependence-aware resolution*

## Central contribution

The method treats scientific attribution as a sequential support-and-resolution problem rather than as predictive validation alone:

1. verify that each declared contrast is estimable from the acquired design;
2. repair missing support at minimum acquisition cost;
3. propagate contrast-aligned numerical uncertainty under an explicitly declared dependence construction; and
4. acquire the next computation until the claim is resolved above the minimum, resolved below the minimum, or retained as indeterminate.

Release `v1.2.0` adds the exhaustive complete-lattice dependence audit. Every archived geometry satisfying the frozen native-cache and mapping-completeness rule was included: four geometries on fine and shared supports, giving eight audits. With identical CFD fields, marginal mapping outputs, response norm, signal functional, and exact finite-set diameter functional, changing only admissible cross-cell coupling reduced the robustness ratio in 8/8 audits and produced seven paired classification reversals at the frozen `R = 3` rule. Independent recombination enlarged the exact finite-set diameter by 10.83-15.71 times. The common diagnostic separation interval is `1.6442596248 < tau < 2.2831676147`.

## Latest reproducibility package

The latest self-contained package is in `v16_method_upgrade/`.

```bash
python -m pip install -r v16_method_upgrade/requirements.txt
cd v16_method_upgrade
python run_all.py
```

The clean-extraction verification runs 15 tests, reconstructs the original four-audit result, and independently rebuilds all eight V16 audits from released float64 weighted-vector caches. The maximum absolute difference from expected output is approximately `2.66e-14`.

## Reproducibility boundary

The released vectors reproduce the numerical decision and dependence-set comparison exactly. They do not provide a lightweight rerun of the underlying OpenFOAM simulations and they do not represent total CFD uncertainty. See `v16_method_upgrade/REPRODUCIBILITY_BOUNDARY.md` before interpretation.

## Repository map

- `v16_method_upgrade/`: current V16 protocols, code, tests, eight float64 audit caches, expected outputs, figures, and hashes.
- `v9_method_upgrade/`: preceding V9 method-upgrade release retained for provenance.
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