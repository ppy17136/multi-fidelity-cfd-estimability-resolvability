"""Closed-loop claim-state acquisition benchmark.

This constructed experiment compares row-oriented design criteria with a policy
that first restores contrast support and then switches to actions that reduce
the declared contrast-level numerical bound.
"""

from __future__ import annotations

import csv
import json
from dataclasses import dataclass
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
SEED = 20260730
REPLICATES = 5000
COEFFICIENTS = np.asarray([1.0, -1.0, -1.0, 1.0])
MINIMUM_EFFECT = 3.0
BUDGET = 30.0


@dataclass(frozen=True)
class Action:
    name: str
    kind: str
    cost: float
    cell: int | None = None
    bound: float | None = None


ACTIONS = [
    Action("acquire_11_base", "row", 4.0, cell=3, bound=0.40),
    Action("acquire_11_refined", "row", 10.0, cell=3, bound=0.10),
    Action("strict_00", "row", 3.0, cell=0, bound=0.10),
    Action("strict_10", "row", 3.5, cell=1, bound=0.10),
    Action("strict_01", "row", 3.5, cell=2, bound=0.10),
    Action("strict_11", "row", 4.0, cell=3, bound=0.10),
    Action("independent_mapping_audit", "evidence", 5.0),
]


def estimable(bounds: np.ndarray) -> bool:
    return bool(np.all(np.isfinite(bounds)))


def estimate_contrast(estimates: np.ndarray) -> float:
    return float(COEFFICIENTS @ estimates)


def contrast_bound(bounds: np.ndarray) -> float:
    if not estimable(bounds):
        return float("inf")
    return float(np.sum(np.abs(COEFFICIENTS) * bounds))


def resolved(estimates: np.ndarray, bounds: np.ndarray) -> bool:
    if not estimable(bounds):
        return False
    lower = max(abs(estimate_contrast(estimates)) - contrast_bound(bounds), 0.0)
    return bool(lower > MINIMUM_EFFECT)


def information_matrix(bounds: np.ndarray) -> np.ndarray:
    info = np.zeros((4, 4))
    for cell, bound in enumerate(bounds):
        if np.isfinite(bound):
            # Uniform(-q,q) has variance q^2/3.
            variance = max(bound * bound / 3.0, 1e-12)
            info[cell, cell] += 1.0 / variance
    return info


def contrast_variance(bounds: np.ndarray) -> float:
    if not estimable(bounds):
        return float("inf")
    return float(COEFFICIENTS @ np.linalg.pinv(information_matrix(bounds)) @ COEFFICIENTS)


def prediction_variance(bounds: np.ndarray) -> float:
    info = information_matrix(bounds)
    ridge = 1e-8
    return float(np.trace(np.linalg.inv(info + ridge * np.eye(4))))


def d_score(bounds: np.ndarray) -> float:
    info = information_matrix(bounds)
    sign, value = np.linalg.slogdet(info + 1e-8 * np.eye(4))
    return float(value if sign > 0 else -np.inf)


def action_available(action: Action, bounds: np.ndarray, used: set[str]) -> bool:
    if action.name in used:
        return False
    if action.name == "strict_11" and not np.isfinite(bounds[3]):
        return False
    if action.kind == "row" and action.cell is not None:
        current = bounds[action.cell]
        if np.isfinite(current) and action.bound is not None and action.bound >= current:
            return False
    if action.kind == "evidence" and not estimable(bounds):
        return False
    return True


def preview_bounds(bounds: np.ndarray, action: Action) -> np.ndarray:
    trial = bounds.copy()
    if action.kind == "row" and action.cell is not None and action.bound is not None:
        trial[action.cell] = min(trial[action.cell], action.bound)
    elif action.name == "independent_mapping_audit":
        trial = np.minimum(trial, 0.20)
    return trial


def choose_action(
    policy: str,
    bounds: np.ndarray,
    used: set[str],
    rng: np.random.Generator,
) -> Action | None:
    available = [a for a in ACTIONS if action_available(a, bounds, used)]
    if not available:
        return None

    if policy in {"proposed_state_switch", "c_optimal_then_resolution"} and estimable(
        bounds
    ):
        current = contrast_bound(bounds)
        scored = []
        for action in available:
            trial = preview_bounds(bounds, action)
            reduction = current - contrast_bound(trial)
            scored.append((-reduction / action.cost, action.cost, action.name, action))
        scored.sort()
        return scored[0][-1]

    if policy == "proposed_state_switch" and not estimable(bounds):
        support_actions = [
            a
            for a in available
            if a.kind == "row" and a.cell == 3 and np.isfinite(preview_bounds(bounds, a)[3])
        ]
        return min(support_actions, key=lambda a: (a.cost, a.name))

    if policy in {"c_optimal", "c_optimal_then_resolution"}:
        scored = []
        for action in available:
            trial = preview_bounds(bounds, action)
            value = contrast_variance(trial)
            objective = value * action.cost
            scored.append((objective, action.cost, action.name, action))
        scored.sort()
        return scored[0][-1]

    if policy == "D_optimal":
        current = d_score(bounds)
        scored = []
        for action in available:
            trial = preview_bounds(bounds, action)
            gain = d_score(trial) - current
            scored.append((-gain / action.cost, action.cost, action.name, action))
        scored.sort()
        return scored[0][-1]

    if policy == "prediction_variance":
        current = prediction_variance(bounds)
        scored = []
        for action in available:
            trial = preview_bounds(bounds, action)
            reduction = current - prediction_variance(trial)
            scored.append((-reduction / action.cost, action.cost, action.name, action))
        scored.sort()
        return scored[0][-1]

    if policy == "high_fidelity_heavy":
        order = [
            "acquire_11_refined",
            "strict_11",
            "strict_10",
            "strict_01",
            "strict_00",
            "independent_mapping_audit",
            "acquire_11_base",
        ]
        lookup = {a.name: a for a in available}
        for name in order:
            if name in lookup:
                return lookup[name]
        return None

    if policy == "random":
        return available[int(rng.integers(0, len(available)))]

    raise ValueError(f"unknown policy {policy}")


def apply_action(
    action: Action,
    true_cells: np.ndarray,
    estimates: np.ndarray,
    bounds: np.ndarray,
    rng: np.random.Generator,
) -> tuple[np.ndarray, np.ndarray]:
    estimates = estimates.copy()
    bounds = bounds.copy()
    if action.kind == "row" and action.cell is not None and action.bound is not None:
        cell = action.cell
        new_bound = action.bound
        estimates[cell] = true_cells[cell] + rng.uniform(-new_bound, new_bound)
        bounds[cell] = min(bounds[cell], new_bound)
    elif action.name == "independent_mapping_audit":
        new_bound = 0.20
        for cell in range(4):
            estimates[cell] = true_cells[cell] + rng.uniform(
                -new_bound, new_bound
            )
        bounds[:] = new_bound
    return estimates, bounds


def run_one(
    policy: str,
    true_delta: float,
    rng: np.random.Generator,
) -> dict:
    true_cells = np.asarray([10.0, 11.0, 9.5, 10.5 + true_delta])
    estimates = np.full(4, np.nan)
    bounds = np.full(4, np.inf)
    for cell in (0, 1, 2):
        bounds[cell] = 0.40
        estimates[cell] = true_cells[cell] + rng.uniform(-0.40, 0.40)

    spent = 0.0
    used: set[str] = set()
    actions: list[str] = []
    while spent < BUDGET and not resolved(estimates, bounds):
        action = choose_action(policy, bounds, used, rng)
        if action is None or spent + action.cost > BUDGET:
            break
        estimates, bounds = apply_action(
            action, true_cells, estimates, bounds, rng
        )
        used.add(action.name)
        actions.append(action.name)
        spent += action.cost

    decision = resolved(estimates, bounds)
    truth = abs(true_delta) > MINIMUM_EFFECT
    return {
        "policy": policy,
        "true_delta": true_delta,
        "truth_above_minimum_effect": truth,
        "resolved": decision,
        "valid_decision": decision == truth,
        "false_resolved": decision and not truth,
        "false_unresolved": (not decision) and truth,
        "cost": spent,
        "contrast_estimate": (
            estimate_contrast(estimates) if estimable(bounds) else float("nan")
        ),
        "contrast_bound": contrast_bound(bounds),
        "action_count": len(actions),
        "actions": ";".join(actions),
    }


def bootstrap_median_ci(values: np.ndarray, rng: np.random.Generator) -> list[float]:
    indices = rng.integers(0, len(values), size=(20_000, len(values)))
    medians = np.median(values[indices], axis=1)
    return [float(x) for x in np.quantile(medians, [0.025, 0.5, 0.975])]


def main() -> None:
    policies = [
        "proposed_state_switch",
        "c_optimal_then_resolution",
        "c_optimal",
        "D_optimal",
        "prediction_variance",
        "high_fidelity_heavy",
        "random",
    ]
    true_deltas = [2.0, 5.0]
    master = np.random.SeedSequence(SEED)
    children = master.spawn(len(policies) * len(true_deltas) * REPLICATES)
    rows = []
    k = 0
    for delta in true_deltas:
        for policy in policies:
            for replicate in range(REPLICATES):
                rng = np.random.default_rng(children[k])
                k += 1
                row = run_one(policy, delta, rng)
                row["replicate"] = replicate
                rows.append(row)

    with (ROOT / "04_closed_loop_acquisition_results.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    bootstrap_rng = np.random.default_rng(SEED + 1)
    for delta in true_deltas:
        for policy in policies:
            subset = [r for r in rows if r["true_delta"] == delta and r["policy"] == policy]
            resolved_rows = [r for r in subset if r["resolved"]]
            costs = np.asarray([r["cost"] for r in resolved_rows], dtype=float)
            summary_rows.append(
                {
                    "true_delta": delta,
                    "policy": policy,
                    "resolved_fraction": float(np.mean([r["resolved"] for r in subset])),
                    "valid_decision_fraction": float(
                        np.mean([r["valid_decision"] for r in subset])
                    ),
                    "false_resolved_fraction": float(
                        np.mean([r["false_resolved"] for r in subset])
                    ),
                    "false_unresolved_fraction": float(
                        np.mean([r["false_unresolved"] for r in subset])
                    ),
                    "median_cost_to_resolved": (
                        float(np.median(costs)) if len(costs) else None
                    ),
                    "median_cost_to_resolved_bootstrap_95": (
                        bootstrap_median_ci(costs, bootstrap_rng)
                        if len(costs)
                        else None
                    ),
                    "most_common_action_path": max(
                        {r["actions"] for r in subset},
                        key=lambda path: sum(r["actions"] == path for r in subset),
                    ),
                }
            )

    summary = {
        "protocol": "V9-20260730",
        "seed": SEED,
        "replicates_per_policy_scenario": REPLICATES,
        "minimum_effect": MINIMUM_EFFECT,
        "budget": BUDGET,
        "rows": summary_rows,
        "claim_boundary": (
            "This is a constructed decision benchmark.  It tests acquisition "
            "logic and error-control behaviour, not external CFD validity."
        ),
    }
    (ROOT / "04_closed_loop_acquisition_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

