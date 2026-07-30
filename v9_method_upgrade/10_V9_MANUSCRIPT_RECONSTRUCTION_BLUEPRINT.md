# V9 manuscript reconstruction blueprint

## Recommended title

**Claim-adaptive multi-fidelity design: Minimum-cost support repair and
uncertainty-set-aware resolution control for CFD**

The title leads with a computational method, not a reporting principle.

## One-sentence method claim

We formulate multi-fidelity attribution as a sequential computational-design
problem that first acquires the minimum-cost support needed to estimate a
predeclared contrast and then switches to actions that reduce contrast-level
numerical uncertainty until a valid effect-present, effect-absent, or
indeterminate decision is reached.

## Editor-visible discovery

**The geometry of the numerical uncertainty set can reverse an attribution
decision even when its marginal cellwise scales are unchanged.**

In the real CFD audit, all four fine/shared-support combinations were robust
across a frozen ensemble of complete mapping pipelines, whereas all four
failed an adversarial independent-cell envelope. The joint mapping diameter
was only 7.2%-8.1% of the axis-aligned envelope because the joint construction
preserved common-mode cancellation across fidelity cells.

## Four contributions

1. A general minimum-cost support-repair optimization for arbitrary linear
   contrasts, multilevel/multi-axis designs, unequal costs, and multiple
   simultaneous targets, with exact and scalable greedy solvers.
2. Deterministic bounded-error and covariance-based contrast-level resolution
   rules, with assumptions and error-rate statements kept separate.
3. A state-dependent acquisition algorithm that switches among support,
   mapping, and repeatability actions according to the current decision
   bottleneck.
4. Real CFD evidence showing both the benefit and the limitation of the
   method: mapping-ensemble robustness depends on uncertainty-set geometry,
   while a second obstruction-flow application transfers Gate 1 but honestly
   leaves Gate 2 unassessed.

## Main-paper structure

1. Introduction
2. Claim-adaptive multi-fidelity design
   1. Target contrast and design-supported estimability
   2. Minimum-cost exact and greedy support repair
   3. Contrast-level deterministic and probabilistic resolution
   4. Joint versus axis-aligned numerical uncertainty sets
   5. State-dependent acquisition and stopping rule
3. Verification and benchmark protocol
   1. Support-repair suites
   2. Gate-2 calibration
   3. Prediction-attribution diagnostic
   4. Nontrivial 3x2 closed-loop benchmark and baselines
4. CFD applications
   1. Periodic-hill field contrast
   2. Frozen complete mapping ensemble
   3. Obstruction-flow scalar transfer audit
5. Results
   1. Repair optimality and scalability
   2. False-decision calibration
   3. Cost to valid decision
   4. Real-CFD uncertainty-set reversal
   5. Transfer and unresolved boundary
6. Discussion
7. Conclusions

## Main figures

1. Algorithm diagram: support state -> resolution state -> stopping decision.
2. Minimum-cost repair benchmark and exact/greedy gap.
3. Cost-to-valid-decision curves for proposed, c-/D-/A-optimal,
   prediction-variance, high-fidelity-heavy, random, and ablation policies.
4. Gate-2 calibration: false-present/false-absent rates and decision regions.
5. Real CFD uncertainty geometry: joint mapping diameter versus axis-aligned
   envelope on fine/shared support.
6. Cross-physics transfer table/compact panel for pressure-loss contrasts.

## Wording rules

- Use `design-supported estimability under the declared model`, not
  `impossible to estimate` without qualification.
- Use `robust within the prespecified mapping ensemble`, not `numerically
  resolved` without scope.
- Use `adversarial independent-cell envelope`, not `total CFD uncertainty`.
- Use `verification` for known-truth experiments.
- Do not claim universal superiority over optimal design.
- Do not present a factor of three as a confidence level.
- Do not retain a section titled `Editor-visible novelty`.

## Placement of legacy results

- Keep the 480-pair prediction-attribution result as motivation, not the main
  contribution.
- Keep the 48-replicate equal-cost result with bootstrap intervals as a
  supporting benchmark.
- Move exhaustive 2x2 enumeration to a unit-test/supplement role.
- Keep 79 verified references and preserve the two required institutional
  citations; renumber only after the V9 body is stable.

