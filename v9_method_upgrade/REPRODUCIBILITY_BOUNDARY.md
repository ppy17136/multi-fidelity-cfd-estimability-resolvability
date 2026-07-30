# V9 reproducibility boundary

This directory freezes the method-upgrade evidence for the claim-adaptive
multi-fidelity design manuscript.

## Directly rerunnable components

The following commands are self-contained within this directory:

```bash
python 01_run_support_repair_benchmarks.py
python 02_run_gate2_calibration.py
python 04b_run_multilevel_closed_loop_acquisition.py
python -m pytest -q
python 08_build_method_upgrade_evidence_summary.py
python build_v9_main_figures.py
```

The seeded simulations regenerate the support-repair, resolution-calibration,
and closed-loop decision metrics. Wall-clock timings, floating-point residuals
below numerical tolerance, and the ordering of equally frequent tied action
paths can vary with BLAS, Python, and hardware; the archived tables are the
frozen manuscript records. The unit-test suite checks support feasibility and
tie-breaking, contrast decisions, multilevel stopping, the frozen joint-mapping
outputs, and the second-physics transfer outputs.

## Archived-output components

Scripts `03`, `05`, `06`, `07`, and `09` document analyses that consume prior
CFD-derived fields or archived project tables. Their frozen CSV/JSON outputs,
protocols, implementation records, hashes, and decision tables are included
here. Raw CFD time histories, processor directories, and third-party DNS
payloads are not redistributed. Consequently, the archived decisions and main
figures are verifiable from this release, whereas end-to-end regeneration of
the raw mapping fields requires the original CFD outputs and public DNS data
identified in the parent repository.

## Claim boundaries

- Greedy repair feasibility is not a claim of cost optimality.
- Normal-model calibration is not transferred to deterministic CFD sensitivity
  runs without a justified stochastic model.
- Robustness across four frozen mapping pipelines is not total CFD uncertainty.
- The pressure-loss interaction is estimable but not declared resolved because
  no independent contrast-level error bound was available.
- Closed-loop cost rankings are benchmark-specific, not universal dominance.
