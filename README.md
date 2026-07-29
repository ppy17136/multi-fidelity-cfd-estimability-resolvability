# Ordered two-gate estimability-resolvability rule in multi-fidelity CFD

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.21670992.svg)](https://doi.org/10.5281/zenodo.21670992)

This repository contains the reproducibility materials for:

> *Estimable but unresolved: completing the fidelity lattice is not enough for
> interaction attribution in multi-fidelity CFD*

## Central result

Completing a multi-fidelity factorial lattice can make an interaction contrast
algebraically estimable without making it numerically resolvable. The study
therefore applies an ordered rule:

1. **Gate 1 — contrast estimability:** verify that the target contrast lies in
   the row space of the acquired design.
2. **Gate 2 — numerical resolvability:** require the observed contrast to
   exceed a prospectively frozen, response-aligned numerical evidence floor.

The manufactured benchmark contains known negative, transition, and positive
controls. In the real-CFD demonstration, all four complete lattices pass Gate
1, while all four interaction-to-floor ratios remain below the frozen Gate-2
threshold of 3.

## Repository map

- `theory/`: row-space derivation and interpretation.
- `protocols/`: prospectively frozen decision and implementation records.
- `code/`: manufactured benchmark, audit, verification, and figure scripts.
- `data/manufactured/`: deterministic seeded benchmark outputs.
- `data/derived/`: compact real-CFD audit results and scientific decisions.
- `data/derived_arrays/`: compact authored CFD-derived arrays.
- `data/manifests/`: portable case and cache manifests without local paths.
- `figures/`: publication figures in PNG and SVG formats.
- `provenance/`: completion markers and SHA-256 inventory.

## Quick verification

From the repository root:

```bash
python code/verify_two_gate_results.py
```

To regenerate the manufactured benchmark and Figure 1:

```bash
python code/199_run_manufactured_two_gate_benchmark.py
```

To regenerate the real-CFD summary Figure 2:

```bash
python code/204_build_real_cfd_two_gate_figure.py
```

The two figure scripts are self-contained apart from the Python dependencies
listed in `requirements.txt`. They write outputs beside the scripts; compare
the regenerated files or values with the archived copies in `figures/`.

## Data boundaries

This archive contains authored code, protocols, derived CFD audit data,
compact derived arrays, figures, hashes, and completion evidence. It does not
contain:

- downloaded journal articles;
- third-party public DNS payloads;
- raw processor directories or full CFD time histories;
- machine-specific runtime files or absolute paths.

Public periodic-hill DNS data must be obtained from and credited to the
original provider:

H. Xiao, J.-L. Wu, S. Laizet, and L. Duan, “Flows over periodic hills of
parameterized geometries: A dataset for data-driven turbulence modeling from
direct simulations,” *Computers & Fluids* 200 (2020) 104431.
https://doi.org/10.1016/j.compfluid.2020.104431.

## Citation and DOI

Repository: https://github.com/ppy17136/multi-fidelity-cfd-estimability-resolvability

- Version v1.0.0 DOI: https://doi.org/10.5281/zenodo.21670993
- Concept DOI for all versions: https://doi.org/10.5281/zenodo.21670992

The version DOI identifies the immutable v1.0.0 archive. The concept DOI
always resolves to the latest archived version.

## Licences

- Source code is licensed under the MIT License; see `LICENSE-CODE`.
- Authored data, figures, and documentation are licensed under the Creative
  Commons Attribution 4.0 International License; see `LICENSE-DATA`.
- Third-party data are not redistributed here and remain subject to the terms
  of their original providers.

## Authorship and funding

The author order follows the associated article.
This work was supported by the Youth Project of the Liaoning Education
Department of China under Grant No. LJKQZ20222277.
