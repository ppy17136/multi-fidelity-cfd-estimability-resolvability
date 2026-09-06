# Auxiliary robust-union analyses

These results evaluate functionals on the union of all pathway assignments
whose disagreement cost does not exceed the budget. This is a different
object from optimizing functionals evaluated inside one common partition.
They are not the source of the main CFD numerical claims or Figures 7.1/7.2.
The manuscript-authoritative CFD profiles and summary are the root 37/42 files.

From the package root run `python run_all.py --include-robust-union`.
The default command skips this optional group. Writers run in an external copy.
Output comparisons are grouped separately in RECONSTRUCTION_REPORT.json.

28 uses unit edge costs; 29 derives diameter-tolerance costs from 28;
33 derives robust-union ratio threshold frontiers from 28; 35 uses normalized
axis-cost scenarios. These budgets need not have the same scale. This ratio
classification is not the scalar decision certificate defined in the paper.
The CSV/JSON payloads are byte-identical to the v2.1.7 auxiliary outputs;
their filenames and locations now identify the robust-union object explicitly.
Builders change only paths, names, and explanatory wording, not computations.
