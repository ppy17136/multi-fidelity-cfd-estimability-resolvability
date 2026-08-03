# Portable eight-audit reconstruction

Run from the repository root:

```text
python V16_complete_lattice_mapping_audit/run_v16_audit.py
```

The command reads only the relative files under this directory and rebuilds
all eight coherent-versus-independent mapping-dependence audits from released
float64 weighted cell vectors.

The `input/complete_lattice_index.csv` file records the frozen completeness
scan. Four geometries meet the rule, all four are included, and each is audited
on the fine F100 and shared F000 supports. No ratio or classification enters the
eligibility rule.

The released vectors are derived data, not raw CFD fields. They preserve all
information required for the finite-set quantities reported in the manuscript:
the minimum weighted norm, exact set diameter, robustness ratio, frozen R=3
classification, and diagnostic common separation interval.
