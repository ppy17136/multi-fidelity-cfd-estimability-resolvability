# Claim-adaptive multi-fidelity design for CFD

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21670992.svg)](https://doi.org/10.5281/zenodo.21670992)

This repository contains the reproducibility materials for:

> *Claim-adaptive multi-fidelity design: Minimum-cost support repair and
> uncertainty-set-aware resolution control for CFD*

## Central contribution

The method treats scientific attribution as a sequential support-and-resolution
problem rather than as predictive validation alone:

1. verify that each declared contrast is estimable from the acquired design;
2. repair missing support at minimum acquisition cost;
3. propagate contrast-aligned numerical uncertainty under an explicitly
   declared dependence model; and
4. acquire the next computation until the claim is effect present, effect
   absent, or indeterminate.

The V9 evidence includes six support-repair benchmarks, bounded and
covariance-based Gate-2 calibration, a 3x2 closed-loop acquisition comparison,
a real-CFD joint-versus-independent mapping audit, and a pressure-loss transfer
audit. The real-CFD finding is deliberately scoped: the same four complete
mapping outputs pass under a coherent joint ensemble and fail under an
independent-cell envelope. This demonstrates dependence-model sensitivity of
the attribution decision; it is not a claim about total CFD uncertainty.

## Repository map

- `v9_method_upgrade/`: frozen protocols, algorithms, tests, result tables, and
  the evidence summary for the current method upgrade.
- `figures/v9_method_upgrade/`: current principal figures in PNG and PDF.
- `theory/`, `protocols/`, `code/`, `data/`, `figures/`, and `provenance/`:
  the preceding ordered two-gate release and its compact CFD-derived evidence.
- `THIRD_PARTY_DATA_NOTICE.md`: data ownership and redistribution boundaries.

## Quick verification

```bash
python -m pip install -r requirements.txt
cd v9_method_upgrade
python -m pytest -q
python 01_run_support_repair_benchmarks.py
python 02_run_gate2_calibration.py
python 04b_run_multilevel_closed_loop_acquisition.py
python 08_build_method_upgrade_evidence_summary.py
python build_v9_main_figures.py
```

The current frozen test suite passes 12/12 tests. See
`v9_method_upgrade/REPRODUCIBILITY_BOUNDARY.md` before interpreting any CFD
mapping or transfer result.

## Data boundaries

This archive contains authored code, protocols, compact derived CFD audit data,
figures, hashes, and completion evidence. It does not contain downloaded
articles, third-party DNS payloads, raw processor directories, or full CFD time
histories. Public periodic-hill DNS data remain with the original provider:

H. Xiao, J.-L. Wu, S. Laizet, and L. Duan, "Flows over periodic hills of
parameterized geometries: A dataset for data-driven turbulence modeling from
direct simulations," *Computers & Fluids* 200 (2020) 104431.
https://doi.org/10.1016/j.compfluid.2020.104431.

## Citation and DOI

Repository: https://github.com/ppy17136/multi-fidelity-cfd-estimability-resolvability

- Concept DOI for all versions: https://doi.org/10.5281/zenodo.21670992
- The immutable version DOI is shown on the corresponding GitHub/Zenodo release.

## Licences

- Source code: MIT (`LICENSE-CODE`).
- Authored data, figures, and documentation: CC BY 4.0 (`LICENSE-DATA`).
- Third-party data are not redistributed.

## Authorship and funding

The author order follows the associated article. This work was supported by the
Youth Project of the Liaoning Education Department of China under Grant No.
LJKQZ20222277.
