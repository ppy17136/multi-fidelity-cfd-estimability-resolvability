# Beyond predictive accuracy in multi-fidelity simulation

Public reproducibility schema: `V20-20260803`

Concept DOI for all public versions: https://doi.org/10.5281/zenodo.21670992

This archive supports the manuscript *Beyond predictive accuracy in
multi-fidelity simulation: Minimum-cost contrast repair and dependence-aware
decisions*. V20 corrects the closed-loop comparison by removing a duplicate
policy alias and assigning all policies the same action-indexed random scenario
within each truth/replicate pair. Seven distinct policies are therefore
reported. The finite-pool support certificate remains limited to Gate-1
support restoration; it is not a certificate of globally minimum total cost
to a terminal decision.

## Main components

1. exact and greedy finite-pool repair of predeclared contrast support;
2. three-state deterministic and scalar covariance-based resolution;
3. a seven-policy closed-loop benchmark with paired action-indexed scenarios;
4. a 2,000-repeat, 11-point truth-grid audit of the same action table;
5. an exhaustive float64 mapping comparison across four frozen-inventory
   periodic-hill geometries and two supports.

For the mapping audit, coherent pipeline assignments are a subset of all
independent cellwise assignments. The direction of ratio contraction is
therefore structural; the empirical findings are its magnitude, the
10.83--15.71-fold exact finite-set diameter expansion, and the resulting
mapping-robustness threshold crossings under the declared audit.

## Quick verification

```bash
python -m pip install -r requirements.txt
python run_all.py
```

Use `python run_all.py --rebuild-synthetic` to regenerate the support,
calibration, corrected closed-loop, and truth-grid outputs before verification.
All archive member paths use POSIX separators and are tested after clean Linux
extraction.

## Canonical files

- `01*`: finite-pool and randomized support repair.
- `02*`: 11-point three-state Gate-2 calibration.
- `04b*`: corrected seven-policy closed-loop benchmark.
- `04c*`: corrected 11-point truth-grid sensitivity.
- `07b*`: original symmetric 4-versus-256 mapping audit.
- `07c*` and `derived_mapping_vectors_v12/`: original float64 reconstruction.
- `09*`: second-flow Gate-1 transfer.
- `V16_complete_lattice_mapping_audit/`: eight-audit reconstruction, frozen
  12-case inventory index, float64 inputs, and expected outputs.
- `figures/`: figures synchronized with the manuscript materials.

The released weighted vectors reproduce the numerical mapping-set decisions.
The eligibility check verifies consistency with the released frozen inventory;
it does not independently rescan omitted raw OpenFOAM caches. The archive is
therefore not a lightweight rerun of the underlying CFD simulations.

Historical protocol records are retained solely for provenance.

## Licences

- Code: MIT (`LICENSE-CODE`).
- Authored data, figures, and documentation: CC BY 4.0 (`LICENSE-DATA`).
- Third-party data are not redistributed.
