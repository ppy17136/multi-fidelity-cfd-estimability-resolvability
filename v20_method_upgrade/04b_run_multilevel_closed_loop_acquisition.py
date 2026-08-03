"""Nontrivial closed-loop benchmark with support and discrepancy actions.

The 3x2 cell-mean design targets the extreme-level mixed contrast
mu00-mu01-mu20+mu21. Two target cells and one irrelevant middle-level cell
are initially absent. Candidate actions can acquire rows, reduce a bounded
mapping discrepancy, or reduce a bounded repeatability component.
"""

from __future__ import annotations

import csv
import gzip
import json
from dataclasses import dataclass
from functools import lru_cache
from itertools import combinations
from pathlib import Path

import numpy as np

from support_repair import exact_minimum_cost_repair


ROOT = Path(__file__).resolve().parent
SEED = 20260730
REPLICATES = 5000
MINIMUM_EFFECT = 3.0
BUDGET = 40.0
DESIGN = np.eye(6)
CONTRAST = np.asarray([1.0, -1.0, 0.0, 0.0, -1.0, 1.0])
TARGET_CELLS = tuple(np.flatnonzero(CONTRAST))


@dataclass(frozen=True)
class Action:
    name: str
    kind: str
    cost: float
    cell: int | None = None
    mapping_bound: float | None = None
    repeatability_bound: float | None = None


ACTIONS = [
    Action("acquire_01_base", "row", 4.0, 1, 0.40, 0.10),
    Action("acquire_11_distractor", "row", 1.0, 3, 0.25, 0.08),
    Action("acquire_20_base", "row", 3.0, 4, 0.50, 0.10),
    Action("acquire_21_base", "row", 5.0, 5, 0.35, 0.10),
    Action("strict_00", "mapping", 5.0, 0, 0.08, None),
    Action("strict_01", "mapping", 8.0, 1, 0.08, None),
    Action("strict_20", "mapping", 8.0, 4, 0.08, None),
    Action("strict_21", "mapping", 10.0, 5, 0.08, None),
    Action("replicate_00", "repeatability", 1.5, 0, None, 0.025),
    Action("replicate_01", "repeatability", 1.5, 1, None, 0.025),
    Action("replicate_20", "repeatability", 1.5, 4, None, 0.025),
    Action("replicate_21", "repeatability", 1.5, 5, None, 0.025),
    Action("independent_mapping_audit", "mapping_audit", 6.0),
]

# The policies are deliberately distinct. ``support_then_c_optimal`` is the
# only c-optimal baseline: it first restores declared support, then optimizes
# nominal contrast variance. The former duplicate alias ``c_optimal`` was
# removed in V20 rather than being reported as a second baseline.
POLICIES = (
    "proposed_decision_targeted_switch",
    "bound_reduction_ablation",
    "support_then_c_optimal",
    "D_optimal",
    "prediction_variance",
    "high_fidelity_heavy",
    "random",
)


@dataclass(frozen=True)
class Scenario:
    """Action-indexed standardised errors shared by all policies in one replicate."""

    measurement_noise: dict[str, tuple[float, float]]
    random_uniforms: tuple[float, ...]


def scenario_measurement_labels() -> tuple[str, ...]:
    labels = ["initial:0", "initial:2"]
    labels.extend(
        action.name
        for action in ACTIONS
        if action.kind in {"row", "mapping", "repeatability"}
    )
    labels.extend(
        f"independent_mapping_audit:{cell}" for cell in TARGET_CELLS
    )
    return tuple(labels)


def make_scenario(rng: np.random.Generator) -> Scenario:
    """Create a policy-invariant tape of standardised measurement errors."""
    return Scenario(
        measurement_noise={
            label: (float(rng.uniform(-1.0, 1.0)), float(rng.uniform(-1.0, 1.0)))
            for label in scenario_measurement_labels()
        },
        random_uniforms=tuple(float(rng.random()) for _ in range(len(ACTIONS))),
    )



def observed_cells(mapping: np.ndarray) -> list[int]:
    return list(np.flatnonzero(np.isfinite(mapping)))


def estimable(mapping: np.ndarray) -> bool:
    return bool(np.all(np.isfinite(mapping[list(TARGET_CELLS)])))


def total_bounds(mapping: np.ndarray, repeatability: np.ndarray) -> np.ndarray:
    return mapping + repeatability


def contrast_bound(mapping: np.ndarray, repeatability: np.ndarray) -> float:
    if not estimable(mapping):
        return float("inf")
    cells = list(TARGET_CELLS)
    return float(
        np.sum(
            np.abs(CONTRAST[cells])
            * total_bounds(mapping, repeatability)[cells]
        )
    )


def contrast_estimate(estimates: np.ndarray) -> float:
    cells = list(TARGET_CELLS)
    return float(CONTRAST[cells] @ estimates[cells])


def decision_state(
    estimates: np.ndarray,
    mapping: np.ndarray,
    repeatability: np.ndarray,
) -> str:
    if not estimable(mapping):
        return "indeterminate"
    estimate = abs(contrast_estimate(estimates))
    bound = contrast_bound(mapping, repeatability)
    if estimate - bound > MINIMUM_EFFECT:
        return "effect_present"
    if estimate + bound < MINIMUM_EFFECT:
        return "effect_absent"
    return "indeterminate"


def action_available(
    action: Action,
    mapping: np.ndarray,
    repeatability: np.ndarray,
    used: set[str],
) -> bool:
    if action.name in used:
        return False
    if action.kind == "mapping_audit":
        return estimable(mapping)
    if action.cell is None:
        return False
    cell = action.cell
    if action.kind == "row":
        return not np.isfinite(mapping[cell])
    if not np.isfinite(mapping[cell]):
        return False
    if action.kind == "mapping":
        return bool(action.mapping_bound < mapping[cell])
    if action.kind == "repeatability":
        return bool(action.repeatability_bound < repeatability[cell])
    return False


def preview(
    mapping: np.ndarray,
    repeatability: np.ndarray,
    action: Action,
) -> tuple[np.ndarray, np.ndarray]:
    q = mapping.copy()
    r = repeatability.copy()
    if action.kind == "row":
        q[action.cell] = action.mapping_bound
        r[action.cell] = action.repeatability_bound
    elif action.kind == "mapping":
        q[action.cell] = min(q[action.cell], action.mapping_bound)
    elif action.kind == "repeatability":
        r[action.cell] = min(r[action.cell], action.repeatability_bound)
    elif action.kind == "mapping_audit":
        q[list(TARGET_CELLS)] = np.minimum(
            q[list(TARGET_CELLS)], 0.15
        )
    return q, r


def nominal_information(
    mapping: np.ndarray, repeatability: np.ndarray
) -> np.ndarray:
    info = np.zeros((6, 6))
    for cell in observed_cells(mapping):
        variance = max(repeatability[cell] ** 2 / 3.0, 1e-12)
        info[cell, cell] = 1.0 / variance
    return info


def nominal_contrast_variance(
    mapping: np.ndarray, repeatability: np.ndarray
) -> float:
    if not estimable(mapping):
        return float("inf")
    inverse = np.linalg.pinv(nominal_information(mapping, repeatability))
    return float(CONTRAST @ inverse @ CONTRAST)


def d_score(mapping: np.ndarray, repeatability: np.ndarray) -> float:
    info = nominal_information(mapping, repeatability)
    sign, value = np.linalg.slogdet(info + 1e-8 * np.eye(6))
    return float(value if sign > 0 else -np.inf)


def prediction_variance(
    mapping: np.ndarray, repeatability: np.ndarray
) -> float:
    info = nominal_information(mapping, repeatability)
    return float(np.trace(np.linalg.inv(info + 1e-8 * np.eye(6))))


@lru_cache(maxsize=None)
def cached_support_plan_names(
    observed: tuple[int, ...],
    available_names: tuple[str, ...],
) -> tuple[str, ...]:
    action_lookup = {action.name: action for action in ACTIONS}
    support_actions = [
        action_lookup[name]
        for name in available_names
        if name in action_lookup
        and action_lookup[name].kind == "row"
        and action_lookup[name].cell in TARGET_CELLS
    ]
    if not support_actions:
        return ()
    candidates = np.asarray([DESIGN[action.cell] for action in support_actions])
    costs = np.asarray([action.cost for action in support_actions])
    result = exact_minimum_cost_repair(
        DESIGN[list(observed)], candidates, costs, CONTRAST
    )
    if not result.feasible:
        return ()
    return tuple(support_actions[index].name for index in result.selected)


def support_plan(
    mapping: np.ndarray,
    available: list[Action],
) -> list[Action]:
    lookup = {action.name: action for action in available}
    names = cached_support_plan_names(
        tuple(observed_cells(mapping)),
        tuple(sorted(lookup)),
    )
    return [lookup[name] for name in names]


def minimum_cost_resolution_plan(
    estimates: np.ndarray,
    mapping: np.ndarray,
    repeatability: np.ndarray,
    available: list[Action],
    remaining_budget: float,
) -> list[Action]:
    estimate = abs(contrast_estimate(estimates))
    target_bound = abs(estimate - MINIMUM_EFFECT)
    current = contrast_bound(mapping, repeatability)
    candidates = []
    for action in available:
        q, r = preview(mapping, repeatability, action)
        if contrast_bound(q, r) < current:
            candidates.append(action)
    best_key = None
    best_subset: tuple[Action, ...] = ()
    for size in range(1, len(candidates) + 1):
        for subset in combinations(candidates, size):
            cost = sum(action.cost for action in subset)
            if cost > remaining_budget:
                continue
            q = mapping.copy()
            r = repeatability.copy()
            for action in subset:
                q, r = preview(q, r, action)
            final_bound = contrast_bound(q, r)
            if final_bound >= target_bound:
                continue
            names = tuple(sorted(action.name for action in subset))
            key = (cost, final_bound, size, names)
            if best_key is None or key < best_key:
                best_key = key
                best_subset = subset
    return list(best_subset)


def choose_action(
    policy: str,
    estimates: np.ndarray,
    mapping: np.ndarray,
    repeatability: np.ndarray,
    used: set[str],
    scenario: Scenario,
    action_step: int,
    remaining_budget: float,
) -> Action | None:
    available = [
        action
        for action in ACTIONS
        if action_available(action, mapping, repeatability, used)
        and action.cost <= remaining_budget
    ]
    if not available:
        return None

    if policy in {
        "proposed_decision_targeted_switch",
        "bound_reduction_ablation",
        "support_then_c_optimal",
    } and not estimable(mapping):
        plan = support_plan(mapping, available)
        return min(plan, key=lambda action: (action.cost, action.name)) if plan else None

    if policy == "proposed_decision_targeted_switch":
        plan = minimum_cost_resolution_plan(
            estimates,
            mapping,
            repeatability,
            available,
            remaining_budget,
        )
        if plan:
            return min(plan, key=lambda action: (action.cost, action.name))

    if policy in {
        "proposed_decision_targeted_switch",
        "bound_reduction_ablation",
    }:
        current = contrast_bound(mapping, repeatability)
        scored = []
        for action in available:
            q, r = preview(mapping, repeatability, action)
            reduction = current - contrast_bound(q, r)
            scored.append((-reduction / action.cost, action.cost, action.name, action))
        return min(scored)[-1]

    if policy == "support_then_c_optimal":
        current = nominal_contrast_variance(mapping, repeatability)
        scored = []
        for action in available:
            q, r = preview(mapping, repeatability, action)
            reduction = current - nominal_contrast_variance(q, r)
            scored.append((-reduction / action.cost, action.cost, action.name, action))
        return min(scored)[-1]

    if policy == "D_optimal":
        current = d_score(mapping, repeatability)
        scored = []
        for action in available:
            q, r = preview(mapping, repeatability, action)
            gain = d_score(q, r) - current
            scored.append((-gain / action.cost, action.cost, action.name, action))
        return min(scored)[-1]

    if policy == "prediction_variance":
        current = prediction_variance(mapping, repeatability)
        scored = []
        for action in available:
            q, r = preview(mapping, repeatability, action)
            reduction = current - prediction_variance(q, r)
            scored.append((-reduction / action.cost, action.cost, action.name, action))
        return min(scored)[-1]

    if policy == "high_fidelity_heavy":
        order = [
            "strict_21", "strict_20", "strict_01", "strict_00",
            "acquire_21_base", "acquire_20_base", "acquire_01_base",
            "independent_mapping_audit",
        ]
        lookup = {action.name: action for action in available}
        for name in order:
            if name in lookup:
                return lookup[name]
        return min(available, key=lambda action: (action.cost, action.name))

    if policy == "random":
        unit = scenario.random_uniforms[action_step]
        return available[min(int(unit * len(available)), len(available) - 1)]

    raise ValueError(policy)


def draw_estimate(
    truth: float,
    mapping_bound: float,
    repeatability_bound: float,
    noise: tuple[float, float],
) -> float:
    return float(truth + noise[0] * mapping_bound + noise[1] * repeatability_bound)


def apply_action(
    action: Action,
    truth: np.ndarray,
    estimates: np.ndarray,
    mapping: np.ndarray,
    repeatability: np.ndarray,
    scenario: Scenario,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    estimates = estimates.copy()
    mapping = mapping.copy()
    repeatability = repeatability.copy()
    if action.kind == "row":
        cell = action.cell
        mapping[cell] = action.mapping_bound
        repeatability[cell] = action.repeatability_bound
        estimates[cell] = draw_estimate(
            truth[cell], mapping[cell], repeatability[cell],
            scenario.measurement_noise[action.name],
        )
    elif action.kind == "mapping":
        cell = action.cell
        mapping[cell] = min(mapping[cell], action.mapping_bound)
        estimates[cell] = draw_estimate(
            truth[cell], mapping[cell], repeatability[cell],
            scenario.measurement_noise[action.name],
        )
    elif action.kind == "repeatability":
        cell = action.cell
        repeatability[cell] = min(repeatability[cell], action.repeatability_bound)
        estimates[cell] = draw_estimate(
            truth[cell], mapping[cell], repeatability[cell],
            scenario.measurement_noise[action.name],
        )
    elif action.kind == "mapping_audit":
        for cell in TARGET_CELLS:
            mapping[cell] = min(mapping[cell], 0.15)
            estimates[cell] = draw_estimate(
                truth[cell], mapping[cell], repeatability[cell],
                scenario.measurement_noise[f"{action.name}:{cell}"],
            )
    return estimates, mapping, repeatability


def run_one(
    policy: str,
    true_delta: float,
    rng: np.random.Generator | None = None,
    *,
    scenario: Scenario | None = None,
) -> dict:
    """Run one policy under a supplied paired scenario or a fresh local seed."""
    if scenario is None:
        if rng is None:
            raise ValueError("Provide either rng or scenario")
        scenario = make_scenario(rng)
    truth = np.asarray([10.0, 9.0, 8.5, 8.0, 8.0, 7.0 + true_delta])
    estimates = np.full(6, np.nan)
    mapping = np.full(6, np.inf)
    repeatability = np.full(6, np.inf)
    for cell, q, r in ((0, 0.45, 0.10), (2, 0.30, 0.08)):
        mapping[cell] = q
        repeatability[cell] = r
        estimates[cell] = draw_estimate(
            truth[cell], q, r, scenario.measurement_noise[f"initial:{cell}"]
        )

    spent = 0.0
    used: set[str] = set()
    actions: list[str] = []
    state = decision_state(estimates, mapping, repeatability)
    while spent < BUDGET and state == "indeterminate":
        action = choose_action(
            policy, estimates, mapping, repeatability, used, scenario,
            len(actions), BUDGET - spent,
        )
        if action is None:
            break
        estimates, mapping, repeatability = apply_action(
            action, truth, estimates, mapping, repeatability, scenario
        )
        spent += action.cost
        used.add(action.name)
        actions.append(action.name)
        state = decision_state(estimates, mapping, repeatability)

    truth_state = (
        "effect_present" if abs(true_delta) > MINIMUM_EFFECT else "effect_absent"
    )
    return {
        "policy": policy,
        "true_delta": true_delta,
        "replicate_truth_state": truth_state,
        "decision_state": state,
        "valid_decision": state == truth_state,
        "false_present": state == "effect_present" and truth_state == "effect_absent",
        "false_absent": state == "effect_absent" and truth_state == "effect_present",
        "indeterminate": state == "indeterminate",
        "cost": spent,
        "contrast_estimate": (
            contrast_estimate(estimates) if estimable(mapping) else float("nan")
        ),
        "contrast_bound": contrast_bound(mapping, repeatability),
        "action_count": len(actions),
        "actions": ";".join(actions),
    }


def bootstrap_median_ci(values: np.ndarray, rng: np.random.Generator) -> list[float]:
    batches = []
    for _ in range(20):
        indices = rng.integers(0, len(values), size=(250, len(values)))
        batches.append(np.median(values[indices], axis=1))
    medians = np.concatenate(batches)
    return [float(value) for value in np.quantile(medians, [0.025, 0.5, 0.975])]


def main() -> None:
    true_deltas = [2.0, 5.0]
    children = np.random.SeedSequence(SEED).spawn(len(true_deltas) * REPLICATES)
    rows = []
    child = 0
    for true_delta in true_deltas:
        for replicate in range(REPLICATES):
            scenario = make_scenario(np.random.default_rng(children[child]))
            child += 1
            for policy in POLICIES:
                row = run_one(policy, true_delta, scenario=scenario)
                row["replicate"] = replicate
                rows.append(row)

    with gzip.open(
        ROOT / "04b_multilevel_closed_loop_results.csv.gz", "wt", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary_rows = []
    bootstrap_rng = np.random.default_rng(SEED + 2)
    for true_delta in true_deltas:
        for policy in POLICIES:
            subset = [
                row for row in rows
                if row["true_delta"] == true_delta and row["policy"] == policy
            ]
            valid_costs = np.asarray(
                [row["cost"] for row in subset if row["valid_decision"]], dtype=float
            )
            paths, counts = np.unique([row["actions"] for row in subset], return_counts=True)
            summary_rows.append(
                {
                    "true_delta": true_delta,
                    "policy": policy,
                    "valid_decision_fraction": float(np.mean([row["valid_decision"] for row in subset])),
                    "false_present_fraction": float(np.mean([row["false_present"] for row in subset])),
                    "false_absent_fraction": float(np.mean([row["false_absent"] for row in subset])),
                    "indeterminate_fraction": float(np.mean([row["indeterminate"] for row in subset])),
                    "median_cost_to_valid_decision": float(np.median(valid_costs)) if len(valid_costs) else None,
                    "median_cost_bootstrap_95": bootstrap_median_ci(valid_costs, bootstrap_rng) if len(valid_costs) else None,
                    "most_common_action_path": str(paths[np.argmax(counts)]),
                }
            )

    summary = {
        "protocol": "V20-20260803",
        "parent_protocol": "V12-20260730",
        "revision_reason": "Removed a duplicate policy and used paired action-indexed scenario tapes.",
        "design": "3x2 cell-mean extreme-level contrast with a distractor cell",
        "seed": SEED,
        "replicates_per_policy_scenario": REPLICATES,
        "paired_action_indexed_scenarios": True,
        "policies": list(POLICIES),
        "minimum_effect": MINIMUM_EFFECT,
        "budget": BUDGET,
        "rows": summary_rows,
        "checks": {
            "all_false_present_zero": all(row["false_present_fraction"] == 0 for row in summary_rows),
            "all_false_absent_zero": all(row["false_absent_fraction"] == 0 for row in summary_rows),
            "policy_names_unique": len(POLICIES) == len(set(POLICIES)),
        },
        "claim_boundary": "Constructed verification of acquisition logic under declared bounded errors; not external CFD validation.",
    }
    (ROOT / "04b_multilevel_closed_loop_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
