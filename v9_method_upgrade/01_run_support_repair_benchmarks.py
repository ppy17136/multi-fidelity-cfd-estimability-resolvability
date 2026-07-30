"""Run the frozen nontrivial support-repair benchmark suite."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from support_repair import (
    binary_factorial_design,
    contrast_for_terms,
    exact_minimum_cost_repair,
    greedy_support_repair,
)


ROOT = Path(__file__).resolve().parent


def binary_case(
    name: str,
    dimensions: int,
    model_terms: list[tuple[int, ...]],
    requested_terms: list[tuple[int, ...]],
    observed_indices: list[int],
    cost_weights: list[float],
    max_condition: float | None = None,
) -> dict:
    design, corners, terms = binary_factorial_design(
        dimensions, included_terms=model_terms
    )
    contrast = contrast_for_terms(terms, requested_terms)
    observed = design[observed_indices]
    candidate_indices = [i for i in range(len(design)) if i not in observed_indices]
    candidates = design[candidate_indices]
    costs = np.asarray(
        [
            1.0
            + sum(weight for value, weight in zip(corners[i], cost_weights) if value > 0)
            for i in candidate_indices
        ]
    )

    exact = exact_minimum_cost_repair(
        observed,
        candidates,
        costs,
        contrast,
        max_condition=max_condition,
    )
    greedy = greedy_support_repair(
        observed,
        candidates,
        costs,
        contrast,
        max_condition=max_condition,
    )

    def map_selected(local_indices: list[int]) -> list[int]:
        return [candidate_indices[i] for i in local_indices]

    exact_dict = exact.to_dict()
    greedy_dict = greedy.to_dict()
    exact_dict["selected_global"] = map_selected(exact.selected)
    greedy_dict["selected_global"] = map_selected(greedy.selected)
    exact_dict["selected_corners"] = [corners[i] for i in exact_dict["selected_global"]]
    greedy_dict["selected_corners"] = [
        corners[i] for i in greedy_dict["selected_global"]
    ]

    gap = (
        (greedy.added_cost - exact.added_cost) / exact.added_cost
        if exact.feasible and exact.added_cost > 0
        else 0.0
    )
    return {
        "case": name,
        "design_family": f"2^{dimensions}",
        "dimensions": dimensions,
        "row_count": len(design),
        "parameter_count": design.shape[1],
        "observed_count": len(observed_indices),
        "candidate_count": len(candidate_indices),
        "requested_terms": [list(x) for x in requested_terms],
        "observed_indices": observed_indices,
        "candidate_global_indices": candidate_indices,
        "candidate_costs": costs.tolist(),
        "exact": exact_dict,
        "greedy": greedy_dict,
        "greedy_cost_gap_fraction": gap,
        "greedy_exact_cost_match": bool(abs(greedy.added_cost - exact.added_cost) < 1e-12),
    }


def multilevel_case() -> dict:
    # Six cell-mean rows for A in {0,1,2}, B in {0,1}.  The target compares
    # the B effect between the extreme A levels:
    # (A2,B1)-(A2,B0)-(A0,B1)+(A0,B0).
    design = np.eye(6)
    contrast = np.asarray([[1.0, -1.0, 0.0, 0.0, -1.0, 1.0]])
    observed_indices = [0, 1, 2, 3]
    candidate_indices = [4, 5]
    costs = np.asarray([4.0, 9.0])
    exact = exact_minimum_cost_repair(
        design[observed_indices], design[candidate_indices], costs, contrast
    )
    greedy = greedy_support_repair(
        design[observed_indices], design[candidate_indices], costs, contrast
    )
    return {
        "case": "3x2_extreme_level_interaction",
        "design_family": "3x2",
        "dimensions": 2,
        "row_count": 6,
        "parameter_count": 6,
        "observed_count": 4,
        "candidate_count": 2,
        "requested_terms": ["(A2-A0)xB"],
        "observed_indices": observed_indices,
        "candidate_global_indices": candidate_indices,
        "candidate_costs": costs.tolist(),
        "exact": exact.to_dict() | {
            "selected_global": [candidate_indices[i] for i in exact.selected]
        },
        "greedy": greedy.to_dict() | {
            "selected_global": [candidate_indices[i] for i in greedy.selected]
        },
        "greedy_cost_gap_fraction": (
            (greedy.added_cost - exact.added_cost) / exact.added_cost
        ),
        "greedy_exact_cost_match": bool(
            abs(greedy.added_cost - exact.added_cost) < 1e-12
        ),
    }


def main() -> None:
    cases = [
        binary_case(
            "binary_2x2_missing_corner",
            dimensions=2,
            model_terms=[(), (0,), (1,), (0, 1)],
            requested_terms=[(0, 1)],
            observed_indices=[0, 1, 2],
            cost_weights=[2.0, 5.0],
        ),
        binary_case(
            "binary_2x2_unequal_cost_alternative_support",
            dimensions=2,
            model_terms=[(), (0,), (1,)],
            requested_terms=[(0,), (1,)],
            observed_indices=[0],
            cost_weights=[2.0, 7.0],
        ),
        binary_case(
            "binary_2x2_condition_constrained",
            dimensions=2,
            model_terms=[(), (0,), (1,)],
            requested_terms=[(0,), (1,)],
            observed_indices=[0],
            cost_weights=[1.0, 1.5],
            max_condition=3.0,
        ),
        binary_case(
            "binary_2x3_fractional_pair_contrasts",
            dimensions=3,
            model_terms=[
                (),
                (0,),
                (1,),
                (2,),
                (0, 1),
                (0, 2),
                (1, 2),
            ],
            requested_terms=[(0, 1), (0, 2)],
            observed_indices=[0, 1, 2, 4],
            cost_weights=[1.0, 2.0, 6.0],
        ),
        binary_case(
            "binary_2x4_sparse_two_way_model",
            dimensions=4,
            model_terms=[
                (),
                (0,),
                (1,),
                (2,),
                (3,),
                (0, 1),
                (0, 2),
                (1, 3),
            ],
            requested_terms=[(0, 1), (0, 2), (1, 3)],
            observed_indices=[0, 1, 2, 4, 8],
            cost_weights=[1.0, 2.0, 4.0, 8.0],
        ),
        multilevel_case(),
    ]

    rows = []
    for case in cases:
        for method in ("exact", "greedy"):
            result = case[method]
            rows.append(
                {
                    "case": case["case"],
                    "design_family": case["design_family"],
                    "method": method,
                    "feasible": result["feasible"],
                    "observed_count": case["observed_count"],
                    "candidate_count": case["candidate_count"],
                    "selected_count": len(result["selected"]),
                    "selected_global": json.dumps(result["selected_global"]),
                    "added_cost": result["added_cost"],
                    "rowspace_residual": result["rowspace_residual"],
                    "contrast_variance": result["contrast_variance"],
                    "condition_number": result["condition_number"],
                    "subsets_evaluated": result["subsets_evaluated"],
                    "runtime_s": result["runtime_s"],
                    "greedy_cost_gap_fraction": (
                        case["greedy_cost_gap_fraction"] if method == "greedy" else 0.0
                    ),
                }
            )

    with (ROOT / "01_support_repair_benchmark_results.json").open(
        "w", encoding="utf-8"
    ) as handle:
        json.dump({"protocol": "V9-20260730", "cases": cases}, handle, indent=2)

    with (ROOT / "01_support_repair_benchmark_results.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    exact_matches = sum(case["greedy_exact_cost_match"] for case in cases)
    summary = {
        "protocol": "V9-20260730",
        "case_count": len(cases),
        "exact_feasible_count": sum(case["exact"]["feasible"] for case in cases),
        "greedy_feasible_count": sum(case["greedy"]["feasible"] for case in cases),
        "greedy_exact_cost_matches": exact_matches,
        "greedy_exact_cost_match_fraction": exact_matches / len(cases),
        "maximum_greedy_cost_gap_fraction": max(
            case["greedy_cost_gap_fraction"] for case in cases
        ),
    }
    (ROOT / "01_support_repair_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

