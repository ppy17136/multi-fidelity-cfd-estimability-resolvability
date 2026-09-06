# Minimum-relaxation certificates for pathway dependence in multifidelity uncertainty quantification

## Current release: v2.1.9

The current verification entry point is `pathway_dependence_v2_1_9/`.
Use the [versioned release](https://github.com/ppy17136/multi-fidelity-cfd-estimability-resolvability/releases/tag/v2.1.9) to identify this correction.

```sh
cd pathway_dependence_v2_1_9
python -m pip install -r requirements.txt
python run_all.py
```

Version 2.1.9 adds exact stored-input guards for scalar crossing witnesses, exact integer-capacity two-label cuts, and an explicit zero-cost state for initially indeterminate scalar cases. The independent diffusion verifier also supports deeply nested Windows work directories. See `RELEASE_V2.1.9.md` for the correction scope.

The verifier checks byte integrity, runs 79 automated tests, reconstructs eight CFD audits, enumerates 128 diffusion-chain partitions, and performs 24 independent diffusion-solver checks. Tests and reconstruction execute in a separate copy without modifying the released inputs. Full generation of the 1,024 diffusion solves is documented separately in the package README.

The common-partition 37/42 outputs supply the manuscript CFD results. Optional 28/29/33/35 robust-union analyses remain in `auxiliary/robust_union/` and are checked with `python run_all.py --include-robust-union`. These constructions are not interchangeable.

## Exact accompanying archive

The accompanying code ZIP is preserved at [`release_archives/pathway_coupling_certificates_v2.1.9_20260907.zip`](release_archives/pathway_coupling_certificates_v2.1.9_20260907.zip).

SHA-256: `9cd8a8478e744a1adbdc21dd11e34e209565ff62fc4867d08c7fe82a79699d9b`

The package manifest covers every payload. Its scientific input arrays and figures are unchanged from version 2.1.8; historical timing records are not a new timing benchmark of this correction.

## Citation and provenance

The preceding [v2.1.8 DOI](https://doi.org/10.5281/zenodo.22540664) identifies the earlier release, not this correction. The [concept DOI](https://doi.org/10.5281/zenodo.21670992) identifies the version series. The earlier [v1.3.0 DOI](https://doi.org/10.5281/zenodo.21767134) identifies the source CFD mapping data.

Earlier directories, releases and source records are retained for provenance. They are not the current verification entry point. `PUBLIC_RELEASE_SHA256.csv` describes the tagged deposit; `CURRENT_CHECKOUT_SHA256.csv` covers the corresponding checkout, excluding the manifest files themselves. Package-specific verification uses `pathway_dependence_v2_1_9/SHA256SUMS.csv`.

## Scope and licences

Graph constraints and weights are declared sensitivity models, not calibrated probabilities or physical error estimates. CFD vectors support finite mapping audits; they do not rerun raw OpenFOAM cases, reconstruct omitted raw-cache eligibility, quantify total CFD uncertainty, or establish physical validation. Floating-point MILP and forest outputs remain distinct from exact-arithmetic optimality proofs.

Code: MIT. Authored data, figures and documentation: CC BY 4.0. Downloaded papers and third-party raw datasets are not redistributed in the current package.
