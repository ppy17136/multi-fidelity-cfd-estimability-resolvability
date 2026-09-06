# Minimum-relaxation certificates for pathway dependence in multifidelity uncertainty quantification

## Scope and contents

Accompanying reproducibility package, version 2.1.8. Contents include graph
certificates, source data, checked witnesses, five manuscript figures, eight
CFD mapping-vector ensembles, and a stochastic diffusion example. The graph
and edge costs define sensitivity models, not calibrated failure probabilities.
The CFD inputs support finite mapping audits, not raw OpenFOAM reruns or
reconstruction of omitted raw-cache eligibility.

## Requirements and verification

Python 3.11 or later is recommended. Install requirements.txt, then run:

    python run_all.py

The manifest is checked first. Tests and writers execute in a separate,
retained temporary copy; its path is printed. Set --work-dir PATH to choose
a new external work location. Keep console logs outside the original archive.
Allow disk space for both the package and reconstruction copy. Runtime depends
on the numerical stack and hardware; measurements are not complexity guarantees.

The command invokes 48 core/weight regression tests, eight theory-extension
tests, six reconstruction-policy tests, eight CFD audits, all 128 diffusion
chain partitions, and 24 independent diffusion-solver checks. It does not
run all 1,024 diffusion solves by default. Success is reported as:

    Pathway-dependence certificates v2.1.8 verification passed

## Main results and optional robust-union analyses

The 37/42 common-partition files supply the main CFD results and figures.
The 28/29/33/35 robust-union outputs are retained under auxiliary/robust_union;
they are not used for the main CFD claims. See docs/RESULT_SOURCE_MAP.md.
Default reconstruction skips them. To include and separately verify them:

    python run_all.py --include-robust-union

Four additional release-routing tests check this separation. With the optional
group enabled, eight auxiliary files are compared separately from the five
manuscript outputs. No cross-object numerical equality is assumed.

## Numerical comparisons

Five manuscript-related regenerated CSV/JSON files use relative tolerance 1e-10 and absolute
tolerance 1e-14. Schema, identities, counts, discrete witnesses,
classifications, and certificate costs must agree. These tolerances do not
change extremizer lists or decision thresholds. Boundary-sensitive discrete
differences cause failure. The original manifest is checked again at the end.
RECONSTRUCTION_REPORT.json is written only in the external work directory.
Downloaded-input byte integrity and regenerated numerical agreement are
different checks. Exact stored-array diameter equality and near-attainment
are separate results, as are common-partition profiles and robust unions.

## Full diffusion generation

    python stochastic_diffusion/run_benchmark.py --work-dir NEW_EXTERNAL_DIRECTORY --full

The launcher verifies and copies original precomputation sources, runs their
pilot/manufactured checks, then executes 1,024 solves. The new directory must
not exist and must be outside this package. Compare numerical arrays; timing
records and container metadata can differ. Singleton extremizer lists make
this a transfer check rather than a difficult multilabel instance.

## Integrity and figure reconstruction

CONTENTS.txt lists files; SHA256SUMS.csv protects all payloads. The unchanged
original precomputation sources are in provenance/diffusion_precomputation.zip,
verified against stochastic_diffusion/results/FREEZE.json. Current PROTOCOL.md
explains the specification and workflow; it is not the original frozen file.
Numerical arrays and benchmark records are unchanged.

Figure builders are make_expansion_figures.py and
plot_juq_provenance_cost_profiles.py; requirements-figures.txt lists additional
dependencies. Execute output-writing scripts in a separate copy. The one-command
verifier arranges this isolation. See docs/REPRODUCIBILITY_BOUNDARY.md for scope.

## Citation and licences

CITATION.cff supplies author metadata. This accompanying version has no assigned
public version DOI. The earlier CFD mapping data are version 1.3.0, DOI
10.5281/zenodo.21767134; that DOI does not identify this new package. The delivery
SHA-256 digest identifies this exact file. Code is MIT-licensed; authored data,
figures, and documentation are CC BY 4.0. See LICENSE-CODE and LICENSE-DATA.
