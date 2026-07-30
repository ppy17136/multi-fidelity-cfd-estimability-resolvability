# Frozen CMAME method-upgrade protocol

Date frozen: 2026-07-30  
Status: prospective for all new algorithmic benchmarks and any new CFD
evidence; retrospective analyses of archived data will be labelled explicitly.

## 1. Objective

Upgrade the existing support-resolution audit into an executable,
claim-targeted multi-fidelity computational-design method. The upgraded method
must:

1. find a minimum-cost acquisition set that makes one or more declared
   contrasts design-supported and estimable;
2. propagate cell-level numerical evidence to the declared contrast;
3. switch acquisition objectives according to whether the active bottleneck is
   support or numerical resolution; and
4. report the cost and error-control properties of the resulting decision.

The existing V8 submission package remains immutable and is not overwritten.

## 2. Claim boundary

The method concerns design-supported, label-level fidelity contrasts. It does
not by itself identify a unique physical cause, eliminate model discrepancy,
or convert a mapping sensitivity into complete numerical uncertainty.

For a linear mean model with design matrix \(X\) and declared contrast matrix
\(C\), support passes only when every row of \(C\) lies in
\(\operatorname{Row}(X)\), under the declared model and rank tolerance.

For a contrast estimator \(\widehat{\theta}=L\mathbf y\), numerical evidence
will be propagated from cell-level perturbations. Deterministic and
probabilistic formulations will be kept separate:

\[
B_C=\sup_{\mathbf e\in\mathcal E}\|L\mathbf e\|_W,
\]

and, only when a covariance model is justified,

\[
\Sigma_C=L\Sigma_eL^\mathsf T.
\]

No probabilistic coverage claim will be made from deterministic sensitivities.

## 3. Frozen algorithmic benchmark suite

The support-repair solver will be tested on:

- binary \(2^2\), \(2^3\), and \(2^4\) factorial designs;
- a \(3\times2\) multi-level design;
- unequal candidate costs;
- one and multiple declared contrasts;
- complete, incomplete, and fractional observed designs.

Two solvers will be implemented:

- exact subset enumeration for finite, moderate candidate pools;
- a scalable greedy rank-residual-per-cost solver with a deterministic
  tie-break based on contrast variance/conditioning.

The following baselines are frozen:

- random acquisition;
- high-fidelity-heavy acquisition where a fidelity ordering exists;
- regularized prediction-variance acquisition;
- D-optimal acquisition;
- c-optimal acquisition;
- proposed claim-state-switching acquisition.

Primary support endpoints:

- feasibility recovery rate;
- exact minimum-cost agreement;
- cost optimality gap;
- final contrast row-space residual;
- condition number;
- runtime.

Primary closed-loop endpoints:

- cost to a valid attribution decision;
- false-resolved rate;
- false-unresolved rate;
- decision coverage;
- final conservative decision margin.

## 4. Frozen Gate-2 verification suite

The manufactured benchmark will use predeclared below-boundary, transition,
and above-boundary contrasts. It will include:

- independent bounded perturbations;
- correlated perturbations;
- heteroscedastic perturbations;
- deliberately underestimated floors by 10%, 25%, and 50%;
- threshold sensitivity over 1.5 to 5.0;
- deterministic-bound and covariance-propagation variants.

Monte Carlo results verify implementation and decision behaviour; they are not
external validation.

## 5. Archived-data analyses

The following analyses use outcomes that already exist and will be labelled
retrospective:

- bootstrap intervals for the 480 matched surrogate-pair audit;
- paired bootstrap and win-probability intervals for the 48 equal-cost
  acquisition replicates;
- threshold sensitivity for the four completed CFD lattices;
- reconstruction of a contrast-level mapping-sensitivity bound from all four
  fidelity cells, if the archived fields permit it.

The prospective frozen primary decisions from Protocols 160, 167, 179, and the
manufactured benchmark remain unchanged.

## 6. Prospective real-CFD repair

The real repair targets the two geometries selected before the original
Protocol-179 outcome:

- `alph15-7929-2024`;
- `alph15-13929-2024`.

Both will be retained to avoid selecting only the geometry with the most
favourable observed ratio.

Candidate evidence actions, in order:

1. reconstruct four-cell contrast-level mapping sensitivity using the archived
   fields and a second independently implemented common-support pathway;
2. quantify solver-continuation/restart sensitivity using stricter convergence
   and an independently initialized or continued solve;
3. if still required, add a systematic refinement action chosen by expected
   reduction of the contrast-level bound per measured CPU-hour.

The primary real-CFD result is not required to switch from unresolved to
resolved. A valid repair may instead demonstrate that the decision remains
unresolved while the evidential bound becomes better specified or narrower.

## 7. Primary manuscript claims permitted after completion

Claims are permitted only if the corresponding checks pass:

1. **Minimum-cost support repair** — exact and greedy solvers verified against
   the frozen benchmark suite.
2. **Contrast-level numerical evidence propagation** — cell-to-contrast
   propagation implemented and its evidential scope stated.
3. **State-dependent acquisition** — proposed policy compared under identical
   candidate pools and costs.
4. **Closed-loop CFD evidence repair** — at least one prospective repair path
   completed and reported regardless of outcome.

The phrases `risk-controlled`, `false-decision-controlled`, and
`numerically resolved` without an evidential qualifier are prohibited unless
their assumptions and coverage are demonstrated.

## 8. Reporting rules

- Manufactured `validation` is renamed `verification`.
- The statistic 3.119 is described as the median within-pair ratio, not as the
  ratio of the two reported medians.
- `Physical interaction` is replaced by `fidelity-source interaction` unless a
  causal statement is independently supported.
- `Cannot estimate` is qualified by the declared design and assumptions.
- Confidence intervals and replicate-level distributions accompany central
  synthetic results.
- Existing school-required references remain only where scientifically
  relevant; citation numbering must remain consecutive.

