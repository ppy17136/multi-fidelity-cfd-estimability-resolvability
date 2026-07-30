# V9 method-upgrade evidence and claim matrix

Protocol: `V9-20260730`

## Editor-critical evidence

| Question | Evidence | Result | Claim boundary |
|---|---|---|---|
| Is support repair a general algorithm? | Six factorial/multilevel, unequal-cost and multi-contrast benchmarks | Exact feasibility 6/6; greedy feasibility 6/6; exact-cost match 5/6; maximum greedy gap 23.8% | Greedy is scalable but not guaranteed optimal |
| Is Gate 2 contrast-level and calibrated? | 200,000 realizations per scenario | Bounded-error rule made no false decisions in its declared model; normal-covariance maximum false-present 2.54% and false-absent 5.29% | Normal-model rates do not apply to deterministic CFD sensitivities |
| Does the acquisition method beat strong baselines? | Frozen 3x2 benchmark, 5,000 replicates per policy and scenario | At true contrast 2, proposed median cost 19.5 versus c-optimal 30.0 (35% lower); at true contrast 5 both cost 13.5 | Frozen benchmark, not universal dominance |
| Was a real Gate-2 mapping repair executed? | Two periodic-hill geometries, fine and shared support, linear plus IDW-4/8/16 | Joint-ensemble pass 4/4, ratios 3.72-9.08; axis-aligned pass 0/4, ratios 0.27-0.69 | Mapping robustness only; not total CFD error |
| Why do the two mapping decisions disagree? | Same field outputs, two uncertainty-set geometries | Joint diameter is only 7.2-8.1% of the axis-aligned envelope | Correlated/common-mode pathway variation is preserved only by the joint ensemble |

## Central methodological finding

**Contrast support and uncertainty-set geometry jointly determine attribution validity.**
Cellwise scalar error floors can reject a contrast even when every prespecified
complete mapping pipeline agrees, because the cellwise envelope destroys the
cross-cell error dependence that cancels in the contrast. Therefore the
uncertainty set—not only its marginal scale—must be declared and audited.

## Go/no-go status for a CMAME-facing V9

| Requirement | Status |
|---|---|
| Method/algorithm rather than only a framework | PASS |
| Optimization beyond the classical row-space condition | PASS |
| Contrast-level deterministic and covariance propagation | PASS |
| Risk/coverage calibration | PASS |
| Real mapping repair | PASS, scope-limited |
| c-optimal and other strong baselines | PASS |
| Cost-to-valid-decision | PASS |
| Nontrivial multilevel design | PASS |
| Second flow class or distinctly different fidelity construction | PASS AS SCOPE-LIMITED TRANSFER AUDIT |
| Editor-decodable Algorithm/Figure 1 | PENDING |

## Second-physics transfer audit

Three obstruction-flow conditions supplied complete medium/fine x
`kEpsilon`/`kOmegaSST` designs for the scalar pressure-loss coefficient.
Their mesh-by-closure contrasts were
0.42%-
0.55%
of the four-cell mean magnitude. Gate 1 was therefore auditable in a distinct
flow and response class. Gate 2 was deliberately not assessed because an
independent contrast-level bound on `K` was not archived for all four cells.

## Permitted wording

- Robust across the four prespecified complete mapping pipelines.
- Unresolved under the independent-cell adversarial envelope.
- The state-switching policy reduced cost in the frozen benchmark.
- The deterministic and probabilistic Gate-2 statements have different assumptions.

## Prohibited wording

- Total CFD uncertainty was quantified.
- The CFD interaction was proven physical or causal.
- The factor of three is a universal confidence threshold.
- The proposed policy universally dominates classical optimal design.
- The greedy support solver is always optimal.
