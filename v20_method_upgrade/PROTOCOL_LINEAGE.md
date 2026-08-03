# Protocol lineage

The public reproducibility schema is `V20-20260803`.

The support, scalar Gate-2, and initial finite-set components remain frozen
under their V12 parent records. V16 added the eight-audit weighted-vector
reconstruction. V20 corrects the synthetic closed-loop comparison by removing
the duplicate `c_optimal` alias and using paired action-indexed random scenario
tapes for all seven distinct policies. V20 does not change the archived CFD
fields, weighted vectors, or mapping-set metrics.

Historical protocol filenames are retained to preserve provenance. The
canonical corrected workflows are
`04b_run_multilevel_closed_loop_acquisition.py` and
`04c_run_truth_grid_closed_loop_audit.py`.
