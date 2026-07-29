# Manufactured two-gate benchmark results

Date: 2026-07-28

## Purpose

This benchmark verifies the logic of the two-gate
estimability-resolvability rule against cases with known truth. It is not a
surrogate for the CFD evidence and does not estimate the real CFD numerical
floor.

## Predeclared construction

The response follows the complete two-factor model

\[
y_{AB}=10+1.5A-0.5B+\Delta AB+\varepsilon_{AB},
\qquad A,B\in\{0,1\}.
\]

Each cell-level numerical perturbation is independently generated within the
deterministic bound

\[
|\varepsilon_{AB}|\le 0.25F.
\]

The four-cell interaction error is therefore bounded by \(F\) through the
triangle inequality. The numerical evidence floor is set to \(F=1\), and the
predeclared gate-2 threshold is

\[
|\widehat\Delta|/F>3.
\]

The random-number seed is 20260728. Each complete-lattice scenario uses
200,000 numerical perturbation realisations.

## Scenarios

### Scenario 1: missing corner with high true signal

- Support: cells \(00,10,01\); cell \(11\) missing.
- True interaction: \(\Delta=5F\).
- Design rank: 3.
- Row-space residual for the interaction coefficient: 1.0.
- Gate 1: fail.
- Gate 2: not evaluated.

This scenario demonstrates that effect magnitude cannot repair structural
non-estimability. A fitted emulator can produce a number at the missing corner,
but the factorial interaction is not identified by the acquired support.

### Scenario 2: complete lattice with below-threshold signal

- Support: complete \(2\times2\) lattice.
- True interaction: \(\Delta=2F\).
- Design rank: 4.
- Row-space residual: \(3.36\times10^{-16}\).
- Gate 1: pass.
- Mean observed ratio: 2.0003.
- Observed ratio range in 200,000 realisations: approximately 1.064-2.956.
- Gate-2 pass rate: 0%.

This is the principal gate-1-pass/gate-2-fail manufactured counterexample.

### Scenario 3: threshold sensitivity

- Support: complete \(2\times2\) lattice.
- True interaction: \(\Delta=3F\).
- Gate 1: pass.
- Mean observed ratio: 3.0006.
- Gate-2 pass rate: 50.1505%.

The threshold case is deliberately not used as a positive or negative control.
It demonstrates that the region near the decision boundary is sensitive to the
realised numerical perturbation and should be reported as a transition region,
not forced into a strong physical conclusion.

### Scenario 4: complete lattice with above-threshold signal

- Support: complete \(2\times2\) lattice.
- True interaction: \(\Delta=5F\).
- Gate 1: pass.
- Mean observed ratio: 5.0007.
- Observed ratio range in 200,000 realisations: approximately 4.033-5.967.
- Gate-2 pass rate: 100%.

This positive control demonstrates that the method does not automatically
rule out all interaction signals.

## Exact bounded-error interpretation

For threshold \(\gamma=3\):

- \(|\Delta|\le2F\) is guaranteed not to pass;
- \(|\Delta|>4F\) is guaranteed to pass;
- \(2F<|\Delta|\le4F\) is a transition region in which the realised numerical
  perturbation can affect the gate-2 decision.

The benchmark Monte Carlo results agree with these deterministic guarantees.

## Relation to the CFD results

The current four completed CFD-lattice ratios are approximately 1.861, 1.976,
2.758, and 2.655, all below the predeclared threshold of 3. They resemble the
manufactured gate-1-pass/gate-2-fail regime, but the manufactured benchmark
does not prove that the CFD perturbations are uniformly bounded or that the
current empirical floor is complete.

The defensible CFD conclusion remains:

> The target interactions are estimable from the completed support but are not
> resolved above the currently audited numerical evidence floor.

It is not:

> The interactions are physically absent.

## Reproducibility artefacts

- `198_estimability_resolvability_theory_20260728.md`
- `199_run_manufactured_two_gate_benchmark.py`
- `199_manufactured_two_gate_summary.csv`
- `199_manufactured_two_gate_quantiles.csv`
- `199_manufactured_two_gate_summary.json`
- `199_fig_manufactured_two_gate_benchmark.png`
- `199_fig_manufactured_two_gate_benchmark.svg`
- `MANUFACTURED_TWO_GATE_BENCHMARK_199_COMPLETE.ok`
