# Protocol lineage

The release schema is `V12-20260730`.

Some experiments were first frozen under parent protocol `V9-20260730` and
were not silently relabelled scientifically. Their copied JSON records carry
both the V12 release schema and a `frozen_parent_protocol` field. Gate-2
three-state calibration and the symmetric 4-versus-256 mapping comparison are
the canonical V12 interpretations. The obsolete V9 binary-decision evidence
summary and the predecessor `04` closed-loop workflow are deliberately absent.

Canonical closed-loop workflow: `04b_run_multilevel_closed_loop_acquisition.py`.
The `04c` workflow is an 11-point sensitivity audit of the same frozen action
table, not a replacement policy or an oracle-regret experiment.
