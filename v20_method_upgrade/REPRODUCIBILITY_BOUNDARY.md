# Reproducibility and data boundary

The archive contains authored algorithms, frozen synthetic benchmarks,
compact derived CFD audit tables, mapping-set outputs, tests, figures, a
frozen 12-case complete-lattice inventory, and eight float64 weighted-vector
caches sufficient to reconstruct the manuscript's central finite-set result.
It excludes downloaded papers, third-party DNS payloads, raw OpenFOAM
processor directories, and full CFD time histories.

The released index records which inventory entries met the frozen completeness
fields, and the audit code checks that every entry marked complete is included.
Because the raw native caches are not redistributed, the lightweight archive
does not independently rebuild those completeness flags from OpenFOAM files.

The periodic-hill result is a finite-set mapping-robustness audit. It is not
total CFD uncertainty, experimental validation, or proof of a physical
mesh-by-closure interaction. The obstructed-flow data establish Gate-1
estimability only because no independent contrast-level Gate-2 bound was
available for all four cells.
