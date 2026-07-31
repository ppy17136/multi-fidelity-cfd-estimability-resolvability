# Contrast-directed multi-fidelity design for CFD

Release schema: `V16-20260731`
Public repository release: `v1.2.0`
Concept DOI for all versions: https://doi.org/10.5281/zenodo.21670992

This archive supports the manuscript *Contrast-directed multi-fidelity design
for CFD: Minimum-cost support repair and dependence-aware resolution*. It
retains the deterministic V12 core and adds a self-contained V16 reconstruction
of every objectively eligible complete-lattice mapping audit.

## Main advances

1. minimum-cost finite-pool repair of predeclared contrast support;
2. three-state deterministic and scalar covariance-based resolution;
3. state-dependent acquisition with a 2,000-repeat 11-point truth-grid audit;
4. an exhaustive float64 real-CFD mapping comparison across four objectively
   eligible geometries and two supports. At `R = 3`, coherent sets pass in 7/8
   audits and independent recombinations in 0/8. For any common diagnostic
   threshold `1.6443 < tau < 2.2832`, all eight audit pairs classify as
   coherent-pass and recombined-fail.

## Quick verification

```bash
python -m pip install -r requirements.txt
python run_all.py
```

`run_all.py` executes the 15 V12 core tests, verifies the original four-audit
float64 reconstruction, and independently rebuilds the V16 eight-audit result
from eight released weighted-vector caches. It checks the frozen eligibility
index, all eight audit pairs, 7/8 and 0/8 pass counts at `R = 3`, directional
ratio reduction in 8/8, the nonempty common separation interval, and numerical
agreement with expected outputs.

## Canonical files

- `01*`: finite-pool and randomized support repair.
- `02*`: 11-point three-state Gate-2 calibration.
- `04b*`: canonical closed-loop benchmark.
- `04c*`: 11-point truth-grid sensitivity.
- `07b*`: original symmetric 4-versus-256 mapping audit.
- `07c*` and `derived_mapping_vectors_v12/`: original float64 reconstruction.
- `09*`: second-flow Gate-1 transfer.
- `V16_complete_lattice_mapping_audit/`: portable exhaustive eight-audit
  reconstruction, frozen inclusion index, float64 inputs, and expected outputs.
- `figures/`: manuscript and supplementary figures.

This V16 complete-lattice extension is included in GitHub release `v1.2.0`.
Zenodo assigns a distinct immutable version DOI to that release while the
concept DOI above resolves to the latest archived version.

## Licences

- Code: MIT (`LICENSE-CODE`).
- Authored data, figures, and documentation: CC BY 4.0 (`LICENSE-DATA`).
- Third-party data are not redistributed.
