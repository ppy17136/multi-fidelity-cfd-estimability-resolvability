# Minimum-relaxation certificates for pathway dependence in multifidelity uncertainty quantification

## Current release: v2.1.8

Public version DOI: [10.5281/zenodo.22540664](https://doi.org/10.5281/zenodo.22540664).

This repository contains versioned reproducibility materials. The current paper and verification entry point are in `pathway_dependence_v2_1_8/`.

```sh
cd pathway_dependence_v2_1_8
python -m pip install -r requirements.txt
python run_all.py
```

The verifier checks byte integrity, runs 66 automated tests, reconstructs eight CFD audits, enumerates 128 diffusion-chain partitions, and performs 24 independent diffusion-solver checks. It operates in an external copy without modifying the released inputs. Full generation of all 1,024 diffusion solves is a separate command documented in the package README.

The common-partition 37/42 outputs supply the main CFD results. Optional 28/29/33/35 robust-union analyses are kept separately under `auxiliary/robust_union/` and checked with `python run_all.py --include-robust-union`. They are different uncertainty constructions and are not substituted for the main results.

## Exact accompanying archive

The journal attachment is preserved byte-for-byte at [`release_archives/pathway_coupling_certificates_v2.1.8_20260906.zip`](release_archives/pathway_coupling_certificates_v2.1.8_20260906.zip).

SHA-256: `f76f5249508cc36d2bde64a053bf6b98811c9f56e19ce331f8278d8b83d50d51`

The nested archive's statements about an unassigned version DOI record its pre-publication build state. This release page and repository-level citation metadata identify the subsequent public archive. Its scientific files and internal manifest are unchanged.

## Citation and provenance

Use [10.5281/zenodo.22540664](https://doi.org/10.5281/zenodo.22540664) for this study's code and certificates. The [concept DOI](https://doi.org/10.5281/zenodo.21670992) identifies the historical version series. The earlier [v1.3.0 DOI](https://doi.org/10.5281/zenodo.21767134) identifies the source CFD mapping data and does not identify the new graph-certificate package.

Earlier directories (`v20_method_upgrade`, `v16_method_upgrade`, `v9_method_upgrade`, and the original top-level code/data material) retain historical results. They are not the current verification entry point. The repository-wide historical manifest predates this release; use `PUBLIC_RELEASE_SHA256.csv` and the current package's `SHA256SUMS.csv` for this release. `PUBLIC_RELEASE_SHA256.csv` describes the immutable tagged deposit; `CURRENT_CHECKOUT_SHA256.csv` covers the current checkout after adding the assigned DOI to these citation documents.

## Scope and licences

Graph constraints and weights are declared sensitivity models, not probabilities or physical error estimates. The CFD vectors support finite mapping audits; they do not rerun raw OpenFOAM cases, reconstruct omitted raw-cache eligibility, quantify total CFD uncertainty, or establish physical validation. Numerical solver witnesses are tolerance-qualified, distinct from symbolic proofs.

Code: MIT. Authored data, figures, and documentation: CC BY 4.0. Downloaded papers and third-party raw datasets are not redistributed in the current package.
