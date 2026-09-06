# Reproducibility boundary

This release supports the finite pathway-coupling results in the associated manuscript.

It verifies graph algorithms, optimization witnesses, derived cost profiles, and the eight weighted-vector CFD mapping audits. It does not reproduce raw OpenFOAM runs, infer raw-cache eligibility, estimate total CFD uncertainty, or validate turbulence physics. The `complete_lattice_index.csv` file is a frozen inventory; `input_provenance.json` records its source metadata boundary.

Exact combinatorial certificates are distinguished from tolerance-qualified floating-point MILP certificates. Exact diameter equality for stored float64 arrays is also distinguished from relative near-attainment at tolerance `1e-7`.