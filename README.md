# Auditing dependence assumptions in multifidelity CFD

## Current verification collection: v3.1.0

Version DOI: <https://doi.org/10.5281/zenodo.23112322>.
All-version concept DOI: <https://doi.org/10.5281/zenodo.21670992>.

The current verification entry point is the downloadable
`MCS_Verification_Records_v3.1.0.zip` asset in the
[v3.1.0 release](https://github.com/ppy17136/multi-fidelity-cfd-estimability-resolvability/releases/tag/v3.1.0),
also archived at the version DOI. Download and extract that asset into a
fresh directory before running:

```text
python verify_all.py
```

The automatically generated GitHub source ZIP is a repository snapshot;
it is not the seven-archive verification collection. Earlier source
directories and releases are retained for provenance, not presented as
the current verification entry point.

The collection contains seven archives and a unified verification index:

- Exact rational same-input conditional-risk comparisons.
- Manufactured Poisson numerical-input stability checks.
- Retrospective four-cell certificates.
- Matched native velocity diagnostics and certificates.
- Fine-grid diagnostics, including failed screens.
- Original six-cell records and continuation checks.
- A small exact-arithmetic certificate-hierarchy check.

The outer SHA-256 manifest covers every payload. The aggregate command
runs the seven saved-record suites in separate fresh temporary directories.
Read the included `README.md` and `Verification_Index_MCS.txt` for dependencies,
record-to-claim mappings, and individual reconstruction boundaries. Python
3.11 or newer, NumPy, and SciPy with `scipy.optimize.milp` support are required
for the complete collection. Numerical-library threads are limited to one.

## Scope and provenance

The records accompany *Auditing dependence assumptions in multifidelity CFD:
Sharp fixed-marginal bounds and minimum-cost identity relaxation*.
They reconstruct conditional numerical-pathway diagnostics, not physical
failure frequencies or certified total-CFD output errors. Identity-release
weights are declared relation weights, not measured CFD acquisition costs.
The workflow does not rerun raw OpenFOAM cases or rebuild omitted raw-cache
eligibility. Eight same-input diagnostic combinations concern one geometry,
not eight independent physical experiments. The manufactured follow-up is
prior-result-informed, not blind validation.

Six archives retain their separately verified bytes. The input-stability
archive carries a documented distribution-only protocol revision and
corresponding verifier metadata; scientific inputs, generating code, and
numerical results are unchanged. Historical titles and run metadata are
retained with explicit reconstruction boundaries.

Earlier version 3.0.1 remains at <https://doi.org/10.5281/zenodo.22866833>;
version 2.1.9 remains at <https://doi.org/10.5281/zenodo.22555851>.
Those immutable records are not overwritten.

## Creator metadata

Current creators: Jia Jun Ma, Bai Lin Lü, Tian Xiang Li, Jia Bao Shang,
and Chao Zhao. Creator metadata was corrected on 2026-10-03. Released
files and their checksums are unchanged; creator lists embedded in the
deposited files reflect the pre-correction build metadata.

## Licences and exclusions

Code: MIT. Authored numerical data and documentation: CC BY 4.0.
The current release contains no manuscript text, submission correspondence,
private reviews, restricted physical-experiment data, third-party articles,
DNS payloads, or complete solver cases.
