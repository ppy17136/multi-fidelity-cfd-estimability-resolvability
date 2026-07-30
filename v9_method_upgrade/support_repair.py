"""Claim-targeted support-repair utilities.

The routines in this module distinguish two operations:

1. recover design-supported estimability of declared linear contrasts; and
2. among support-feasible designs, improve contrast variance/conditioning.

The exact solver enumerates finite candidate subsets and is intended as a
reference implementation for moderate pools.  The greedy solver is a scalable
deterministic approximation.
"""

from __future__ import annotations

from dataclasses import dataclass, asdict
from itertools import combinations
from time import perf_counter
from typing import Iterable, Sequence

import numpy as np


Array = np.ndarray


@dataclass
class RepairResult:
    method: str
    feasible: bool
    selected: list[int]
    added_cost: float
    rowspace_residual: float
    contrast_variance: float
    condition_number: float
    subsets_evaluated: int
    runtime_s: float

    def to_dict(self) -> dict:
        return asdict(self)


def as_contrast_matrix(contrast: Array | Sequence[float]) -> Array:
    c = np.asarray(contrast, dtype=float)
    if c.ndim == 1:
        c = c[None, :]
    if c.ndim != 2:
        raise ValueError("contrast must be a vector or matrix")
    return c


def rowspace_residual(design: Array, contrast: Array) -> float:
    """Maximum relative residual after projection onto Row(design)."""
    x = np.asarray(design, dtype=float)
    c = as_contrast_matrix(contrast)
    if x.ndim != 2 or x.shape[1] != c.shape[1]:
        raise ValueError("design and contrast dimensions are incompatible")
    projector = np.linalg.pinv(x) @ x
    residuals = []
    for row in c:
        denom = max(float(np.linalg.norm(row)), 1.0)
        residuals.append(float(np.linalg.norm(row - row @ projector) / denom))
    return max(residuals, default=0.0)


def is_estimable(design: Array, contrast: Array, tolerance: float = 1e-10) -> bool:
    return rowspace_residual(design, contrast) <= tolerance


def effective_condition_number(design: Array, tolerance: float = 1e-12) -> float:
    """Condition number on the nonzero singular subspace."""
    s = np.linalg.svd(np.asarray(design, dtype=float), compute_uv=False)
    if len(s) == 0 or s[0] == 0:
        return float("inf")
    keep = s > tolerance * s[0]
    if not np.any(keep):
        return float("inf")
    return float(s[keep][0] / s[keep][-1])


def contrast_variance(design: Array, contrast: Array) -> float:
    """Unit-noise generalized least-squares contrast variance trace.

    Returns infinity if the contrast is not estimable.  For an estimable
    contrast this generalized-inverse expression is invariant to the chosen
    coefficient representation.
    """
    x = np.asarray(design, dtype=float)
    c = as_contrast_matrix(contrast)
    if not is_estimable(x, c):
        return float("inf")
    information_pinv = np.linalg.pinv(x.T @ x)
    return float(np.trace(c @ information_pinv @ c.T))


def _assembled_design(
    observed: Array, candidate_rows: Array, selected: Iterable[int]
) -> Array:
    selected = list(selected)
    if not selected:
        return np.asarray(observed, dtype=float)
    return np.vstack([observed, candidate_rows[selected]])


def exact_minimum_cost_repair(
    observed: Array,
    candidate_rows: Array,
    candidate_costs: Sequence[float],
    contrast: Array,
    *,
    tolerance: float = 1e-10,
    max_condition: float | None = None,
) -> RepairResult:
    """Enumerate all candidate subsets and return the lexicographic optimum.

    Objective order:
    1. minimum added cost;
    2. minimum contrast variance;
    3. minimum effective condition number;
    4. minimum number of added rows;
    5. lexicographically smallest candidate-index tuple.
    """
    start = perf_counter()
    observed = np.asarray(observed, dtype=float)
    candidates = np.asarray(candidate_rows, dtype=float)
    costs = np.asarray(candidate_costs, dtype=float)
    c = as_contrast_matrix(contrast)
    if len(candidates) != len(costs):
        raise ValueError("candidate rows and costs must have the same length")
    if np.any(costs < 0):
        raise ValueError("candidate costs must be nonnegative")

    best_key = None
    best_selected: tuple[int, ...] | None = None
    best_metrics = None
    evaluated = 0

    for size in range(len(candidates) + 1):
        for subset in combinations(range(len(candidates)), size):
            evaluated += 1
            x = _assembled_design(observed, candidates, subset)
            residual = rowspace_residual(x, c)
            if residual > tolerance:
                continue
            condition = effective_condition_number(x)
            if max_condition is not None and condition > max_condition:
                continue
            variance = contrast_variance(x, c)
            cost = float(costs[list(subset)].sum()) if subset else 0.0
            key = (cost, variance, condition, size, subset)
            if best_key is None or key < best_key:
                best_key = key
                best_selected = subset
                best_metrics = (residual, variance, condition, cost)

    runtime = perf_counter() - start
    if best_selected is None or best_metrics is None:
        residual = rowspace_residual(observed, c)
        return RepairResult(
            method="exact",
            feasible=False,
            selected=[],
            added_cost=float("inf"),
            rowspace_residual=residual,
            contrast_variance=float("inf"),
            condition_number=effective_condition_number(observed),
            subsets_evaluated=evaluated,
            runtime_s=runtime,
        )

    residual, variance, condition, cost = best_metrics
    return RepairResult(
        method="exact",
        feasible=True,
        selected=list(best_selected),
        added_cost=cost,
        rowspace_residual=residual,
        contrast_variance=variance,
        condition_number=condition,
        subsets_evaluated=evaluated,
        runtime_s=runtime,
    )


def greedy_support_repair(
    observed: Array,
    candidate_rows: Array,
    candidate_costs: Sequence[float],
    contrast: Array,
    *,
    tolerance: float = 1e-10,
    max_condition: float | None = None,
) -> RepairResult:
    """Greedy residual-reduction-per-cost repair with deterministic tie-breaks."""
    start = perf_counter()
    observed = np.asarray(observed, dtype=float)
    candidates = np.asarray(candidate_rows, dtype=float)
    costs = np.asarray(candidate_costs, dtype=float)
    c = as_contrast_matrix(contrast)
    remaining = list(range(len(candidates)))
    selected: list[int] = []
    evaluated = 0

    while True:
        x = _assembled_design(observed, candidates, selected)
        residual = rowspace_residual(x, c)
        condition = effective_condition_number(x)
        support_ok = residual <= tolerance
        condition_ok = max_condition is None or condition <= max_condition
        if support_ok and condition_ok:
            break
        if not remaining:
            return RepairResult(
                method="greedy",
                feasible=False,
                selected=selected,
                added_cost=float(costs[selected].sum()) if selected else 0.0,
                rowspace_residual=residual,
                contrast_variance=contrast_variance(x, c),
                condition_number=condition,
                subsets_evaluated=evaluated,
                runtime_s=perf_counter() - start,
            )

        candidates_scored = []
        for index in remaining:
            evaluated += 1
            trial_selected = selected + [index]
            trial = _assembled_design(observed, candidates, trial_selected)
            trial_residual = rowspace_residual(trial, c)
            residual_gain = max(0.0, residual - trial_residual)
            trial_variance = contrast_variance(trial, c)
            trial_condition = effective_condition_number(trial)
            cost = max(float(costs[index]), np.finfo(float).eps)

            # Primary score repairs support.  Once support is achieved, the
            # secondary terms prefer stable, low-variance representations.
            gain_per_cost = residual_gain / cost
            finite_variance = (
                trial_variance if np.isfinite(trial_variance) else 1e300
            )
            candidates_scored.append(
                (
                    -gain_per_cost,
                    trial_residual,
                    finite_variance,
                    trial_condition,
                    costs[index],
                    index,
                )
            )

        candidates_scored.sort()
        chosen = int(candidates_scored[0][-1])
        selected.append(chosen)
        remaining.remove(chosen)

    x = _assembled_design(observed, candidates, selected)
    return RepairResult(
        method="greedy",
        feasible=True,
        selected=selected,
        added_cost=float(costs[selected].sum()) if selected else 0.0,
        rowspace_residual=rowspace_residual(x, c),
        contrast_variance=contrast_variance(x, c),
        condition_number=effective_condition_number(x),
        subsets_evaluated=evaluated,
        runtime_s=perf_counter() - start,
    )


def binary_factorial_design(
    dimensions: int,
    *,
    included_terms: Sequence[tuple[int, ...]] | None = None,
) -> tuple[Array, list[tuple[int, ...]], list[tuple[int, ...]]]:
    """Return effect-coded binary factorial rows and term labels.

    Corners are tuples in ``{-1, +1}^dimensions``.  Terms are represented by
    tuples of factor indices; the empty tuple denotes the intercept.
    """
    if dimensions < 1:
        raise ValueError("dimensions must be positive")
    corners = [
        tuple(-1 if (mask >> j) & 1 == 0 else 1 for j in range(dimensions))
        for mask in range(2**dimensions)
    ]
    if included_terms is None:
        terms: list[tuple[int, ...]] = [()]
        for order in range(1, dimensions + 1):
            terms.extend(combinations(range(dimensions), order))
    else:
        terms = [tuple(term) for term in included_terms]
        if () not in terms:
            terms = [()] + terms
    rows = []
    for corner in corners:
        row = []
        for term in terms:
            value = 1.0
            for index in term:
                value *= corner[index]
            row.append(value)
        rows.append(row)
    return np.asarray(rows, dtype=float), corners, terms


def contrast_for_terms(
    terms: Sequence[tuple[int, ...]],
    requested_terms: Sequence[tuple[int, ...]],
) -> Array:
    rows = []
    for requested in requested_terms:
        row = np.zeros(len(terms), dtype=float)
        try:
            row[terms.index(tuple(requested))] = 1.0
        except ValueError as exc:
            raise ValueError(f"requested term {requested} is not in the model") from exc
        rows.append(row)
    return np.asarray(rows)

