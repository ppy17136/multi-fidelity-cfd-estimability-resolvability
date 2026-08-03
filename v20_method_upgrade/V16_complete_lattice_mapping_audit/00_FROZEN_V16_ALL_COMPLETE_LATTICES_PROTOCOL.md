# Frozen V16 complete-lattice mapping-dependence audit

Eligibility was fixed before inspection of any robustness ratio. A geometry is
eligible only when F000, F010, F100, and F110 have an indexed native CFD cache
and non-empty C, U, and Vc fields. Every eligible geometry is retained.

Each eligible geometry is evaluated on the fine F100 and shared F000 supports.
The four frozen pathways are periodic linear interpolation and periodic IDW
with 4, 8, and 16 neighbours. Interaction coefficients are (+1,-1,-1,+1) for
(F000,F100,F010,F110). The coherent set applies one pathway to all four cells;
the independent set contains all 4^4 cellwise recombinations.

Both sets use S(U)=min ||C||_W, D(U)=max ||C-C'||_W, and R(U)=S(U)/D(U),
computed from float64 weighted vectors. R>=3 is the frozen display rule. The
common interval reported in the manuscript is a diagnostic separation interval,
not a replacement threshold selected after inspecting the data.
