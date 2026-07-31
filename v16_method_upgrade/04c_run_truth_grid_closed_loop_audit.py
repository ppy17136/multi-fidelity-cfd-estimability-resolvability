"""Truth-grid audit of the frozen multilevel closed-loop action table."""

from __future__ import annotations

import csv
import gzip
import importlib.util
import json
import sys
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
BASE_PATH = ROOT / "04b_run_multilevel_closed_loop_acquisition.py"
SEED = 20260731
REPLICATES = 2000
TRUE_DELTAS = (0.0, 1.0, 2.0, 2.5, 2.9, 3.0, 3.1, 3.5, 4.0, 5.0, 6.0)
POLICIES = (
    "proposed_decision_targeted_switch",
    "bound_reduction_ablation",
    "support_then_c_optimal",
    "c_optimal",
    "D_optimal",
    "prediction_variance",
    "high_fidelity_heavy",
    "random",
)


def load_base():
    spec = importlib.util.spec_from_file_location("closed_loop_v10", BASE_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def expected_state(delta: float, threshold: float) -> str:
    magnitude = abs(delta)
    if magnitude > threshold:
        return "effect_present"
    if magnitude < threshold:
        return "effect_absent"
    return "indeterminate"


def bootstrap_median_ci(values, rng):
    if not len(values):
        return None
    batches = []
    for _ in range(10):
        indices = rng.integers(0, len(values), size=(200, len(values)))
        batches.append(np.median(values[indices], axis=1))
    medians = np.concatenate(batches)
    return [float(x) for x in np.quantile(medians, [0.025, 0.5, 0.975])]


def main() -> None:
    base = load_base()
    children = np.random.SeedSequence(SEED).spawn(
        len(TRUE_DELTAS) * len(POLICIES) * REPLICATES
    )
    raw_path = ROOT / "04c_truth_grid_closed_loop_results.csv.gz"
    fieldnames = [
        "true_delta", "policy", "replicate", "expected_state",
        "decision_state", "correct_outcome", "wrong_decisive",
        "indeterminate", "cost", "contrast_estimate", "contrast_bound",
        "action_count", "actions",
    ]
    rows_by_key = {
        (delta, policy): [] for delta in TRUE_DELTAS for policy in POLICIES
    }
    child = 0
    with gzip.open(raw_path, "wt", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for delta in TRUE_DELTAS:
            target = expected_state(delta, base.MINIMUM_EFFECT)
            for policy in POLICIES:
                for replicate in range(REPLICATES):
                    rng = np.random.default_rng(children[child])
                    child += 1
                    result = base.run_one(policy, delta, rng)
                    state = result["decision_state"]
                    row = {
                        "true_delta": delta,
                        "policy": policy,
                        "replicate": replicate,
                        "expected_state": target,
                        "decision_state": state,
                        "correct_outcome": state == target,
                        "wrong_decisive": (
                            state != "indeterminate" and state != target
                        ),
                        "indeterminate": state == "indeterminate",
                        "cost": result["cost"],
                        "contrast_estimate": result["contrast_estimate"],
                        "contrast_bound": result["contrast_bound"],
                        "action_count": result["action_count"],
                        "actions": result["actions"],
                    }
                    writer.writerow(row)
                    rows_by_key[(delta, policy)].append(row)

    rng = np.random.default_rng(SEED + 1)
    summaries = []
    for delta in TRUE_DELTAS:
        for policy in POLICIES:
            subset = rows_by_key[(delta, policy)]
            correct_costs = np.asarray(
                [x["cost"] for x in subset if x["correct_outcome"]], dtype=float
            )
            all_costs = np.asarray([x["cost"] for x in subset], dtype=float)
            summaries.append({
                "true_delta": delta,
                "policy": policy,
                "expected_state": expected_state(delta, base.MINIMUM_EFFECT),
                "correct_outcome_fraction": float(np.mean(
                    [x["correct_outcome"] for x in subset]
                )),
                "wrong_decisive_fraction": float(np.mean(
                    [x["wrong_decisive"] for x in subset]
                )),
                "indeterminate_fraction": float(np.mean(
                    [x["indeterminate"] for x in subset]
                )),
                "median_cost_to_correct_outcome": (
                    float(np.median(correct_costs)) if len(correct_costs) else None
                ),
                "median_cost_bootstrap_95": bootstrap_median_ci(correct_costs, rng),
                "mean_terminal_cost": float(np.mean(all_costs)),
            })

    with (ROOT / "04c_truth_grid_closed_loop_summary.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
        writer.writeheader()
        writer.writerows(summaries)

    by_key = {(x["true_delta"], x["policy"]): x for x in summaries}
    comparisons = []
    for delta in TRUE_DELTAS:
        proposed = by_key[(delta, POLICIES[0])]
        for baseline in POLICIES[1:]:
            other = by_key[(delta, baseline)]
            pf = proposed["correct_outcome_fraction"]
            bf = other["correct_outcome_fraction"]
            if not np.isclose(pf, bf):
                result = "win" if pf > bf else "loss"
                basis = "correct_outcome_fraction"
            else:
                pc = proposed["median_cost_to_correct_outcome"]
                bc = other["median_cost_to_correct_outcome"]
                if pc is None or bc is None or np.isclose(pc, bc):
                    result = "tie"
                else:
                    result = "win" if pc < bc else "loss"
                basis = "median_cost_given_equal_correct_fraction"
            comparisons.append({
                "true_delta": delta,
                "baseline": baseline,
                "result_for_proposed": result,
                "comparison_basis": basis,
                "proposed_correct_fraction": pf,
                "baseline_correct_fraction": bf,
                "proposed_median_cost": proposed["median_cost_to_correct_outcome"],
                "baseline_median_cost": other["median_cost_to_correct_outcome"],
            })

    payload = {
        "protocol": "V11-20260730",
        "parent_frozen_action_table": "V9-20260730 / 04b",
        "seed": SEED,
        "replicates_per_policy_truth": REPLICATES,
        "true_deltas": list(TRUE_DELTAS),
        "policies": list(POLICIES),
        "minimum_effect": base.MINIMUM_EFFECT,
        "budget": base.BUDGET,
        "summary": summaries,
        "pairwise_comparisons_against_proposed": comparisons,
        "comparison_counts": {
            label: sum(x["result_for_proposed"] == label for x in comparisons)
            for label in ("win", "tie", "loss")
        },
        "all_wrong_decisive_fractions_zero": all(
            x["wrong_decisive_fraction"] == 0 for x in summaries
        ),
        "claim_boundary": (
            "Truth-grid robustness audit under one frozen synthetic action "
            "table; not a random-action benchmark, oracle-regret study, or "
            "external CFD validation."
        ),
    }
    (ROOT / "04c_truth_grid_closed_loop_summary.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
