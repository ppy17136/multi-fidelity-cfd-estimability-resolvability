"""Randomized finite-pool and runtime-scaling benchmarks for support repair.

The randomized comparison uses effect-coded binary factorial rows and models
containing the intercept, all main effects, and all two-factor interactions.
Only instances for which the observed design does not support the declared
contrast but the finite candidate pool does are retained.

The separate runtime experiment is deliberately labelled computational: it
uses generic standardized rows to measure how the greedy implementation scales
with candidate-pool size.  It is not used to make statistical optimality
claims.
"""

from __future__ import annotations

import csv
import json
from itertools import combinations
from pathlib import Path
from time import perf_counter

import numpy as np

from support_repair import (
    binary_factorial_design,
    contrast_for_terms,
    exact_minimum_cost_repair,
    greedy_support_repair,
    is_estimable,
)


ROOT = Path(__file__).resolve().parent
SEED = 20260730
TARGET_INSTANCES = 300


def percentile(values: list[float], q: float) -> float:
    return float(np.percentile(np.asarray(values, dtype=float), q))


def make_random_instance(rng: np.random.Generator, instance_id: int) -> dict:
    for _ in range(20_000):
        dimensions = int(rng.integers(3, 7))
        terms = [()]
        terms.extend((j,) for j in range(dimensions))
        terms.extend(combinations(range(dimensions), 2))
        design, corners, terms = binary_factorial_design(
            dimensions, included_terms=terms
        )

        requested_count = int(rng.integers(1, min(5, dimensions) + 1))
        pair_terms = [term for term in terms if len(term) == 2]
        requested = [
            pair_terms[j]
            for j in rng.choice(
                len(pair_terms), size=min(requested_count, len(pair_terms)), replace=False
            )
        ]
        contrast = contrast_for_terms(terms, requested)

        parameter_count = design.shape[1]
        observed_count = int(
            rng.integers(max(3, parameter_count - 5), parameter_count)
        )
        observed_count = min(observed_count, len(design) - 4)
        permutation = rng.permutation(len(design))
        observed_indices = permutation[:observed_count].tolist()
        remaining = permutation[observed_count:].tolist()
        candidate_count = int(rng.integers(4, min(13, len(remaining)) + 1))
        candidate_indices = remaining[:candidate_count]

        observed = design[observed_indices]
        candidates = design[candidate_indices]
        if is_estimable(observed, contrast):
            continue
        if not is_estimable(np.vstack([observed, candidates]), contrast):
            continue

        cost_sigma = float(rng.choice([0.15, 0.50, 1.00]))
        costs = np.exp(rng.normal(0.0, cost_sigma, size=candidate_count))
        costs = np.round(costs, 4)

        exact = exact_minimum_cost_repair(
            observed, candidates, costs, contrast
        )
        greedy = greedy_support_repair(
            observed, candidates, costs, contrast
        )
        if not exact.feasible:
            continue

        gap = (
            (greedy.added_cost - exact.added_cost) / exact.added_cost
            if greedy.feasible and exact.added_cost > 0
            else float("inf")
        )
        return {
            "instance_id": instance_id,
            "dimensions": dimensions,
            "parameter_count": parameter_count,
            "observed_count": observed_count,
            "candidate_count": candidate_count,
            "requested_contrasts": len(requested),
            "cost_sigma": cost_sigma,
            "exact_feasible": exact.feasible,
            "greedy_feasible": greedy.feasible,
            "exact_cost": exact.added_cost,
            "greedy_cost": greedy.added_cost,
            "greedy_cost_gap_fraction": gap,
            "exact_runtime_s": exact.runtime_s,
            "greedy_runtime_s": greedy.runtime_s,
            "exact_subsets_evaluated": exact.subsets_evaluated,
            "greedy_candidates_evaluated": greedy.subsets_evaluated,
            "exact_selected_count": len(exact.selected),
            "greedy_selected_count": len(greedy.selected),
            "observed_indices": observed_indices,
            "candidate_indices": candidate_indices,
            "requested_terms": [list(term) for term in requested],
            "candidate_costs": costs.tolist(),
            "exact_selected_local": exact.selected,
            "greedy_selected_local": greedy.selected,
            "corners": [list(corner) for corner in corners],
        }
    raise RuntimeError("failed to construct a feasible randomized instance")


def run_randomized_suite(rng: np.random.Generator) -> tuple[list[dict], dict]:
    rows = [make_random_instance(rng, i + 1) for i in range(TARGET_INSTANCES)]
    gaps = [
        row["greedy_cost_gap_fraction"]
        for row in rows
        if np.isfinite(row["greedy_cost_gap_fraction"])
    ]
    exact_times = [row["exact_runtime_s"] for row in rows]
    greedy_times = [row["greedy_runtime_s"] for row in rows]
    summary = {
        "seed": SEED,
        "retained_instances": len(rows),
        "generation_rule": (
            "effect-coded binary factorial rows; intercept, main effects, and "
            "all two-factor interactions; retain Gate-1-failing observed "
            "designs repaired by finite pools of 4-12 candidates"
        ),
        "exact_feasible_fraction": float(
            np.mean([row["exact_feasible"] for row in rows])
        ),
        "greedy_feasible_fraction": float(
            np.mean([row["greedy_feasible"] for row in rows])
        ),
        "greedy_exact_cost_match_fraction": float(
            np.mean([abs(gap) <= 1e-12 for gap in gaps])
        ),
        "greedy_cost_gap_median_fraction": percentile(gaps, 50),
        "greedy_cost_gap_p95_fraction": percentile(gaps, 95),
        "greedy_cost_gap_max_fraction": max(gaps),
        "exact_runtime_median_s": percentile(exact_times, 50),
        "exact_runtime_p95_s": percentile(exact_times, 95),
        "greedy_runtime_median_s": percentile(greedy_times, 50),
        "greedy_runtime_p95_s": percentile(greedy_times, 95),
    }
    return rows, summary


def run_runtime_scaling(rng: np.random.Generator) -> tuple[list[dict], dict]:
    pool_sizes = [16, 32, 64, 128, 256, 512, 1024]
    repeats = 20
    rows: list[dict] = []
    parameter_count = 18
    observed_count = parameter_count - 2
    contrast = np.zeros((1, parameter_count))
    contrast[0, -1] = 1.0

    for pool_size in pool_sizes:
        for repeat in range(repeats):
            observed = rng.normal(size=(observed_count, parameter_count))
            observed /= np.linalg.norm(observed, axis=1, keepdims=True)
            candidates = rng.normal(size=(pool_size, parameter_count))
            candidates /= np.linalg.norm(candidates, axis=1, keepdims=True)
            costs = np.exp(rng.normal(0.0, 0.5, size=pool_size))
            start = perf_counter()
            result = greedy_support_repair(
                observed, candidates, costs, contrast
            )
            elapsed = perf_counter() - start
            rows.append(
                {
                    "pool_size": pool_size,
                    "repeat": repeat + 1,
                    "feasible": result.feasible,
                    "selected_count": len(result.selected),
                    "candidates_evaluated": result.subsets_evaluated,
                    "runtime_s": elapsed,
                }
            )

    grouped = {}
    for pool_size in pool_sizes:
        subset = [row for row in rows if row["pool_size"] == pool_size]
        times = [row["runtime_s"] for row in subset]
        grouped[str(pool_size)] = {
            "repeats": len(subset),
            "feasible_fraction": float(np.mean([row["feasible"] for row in subset])),
            "runtime_median_s": percentile(times, 50),
            "runtime_p95_s": percentile(times, 95),
            "selected_count_median": percentile(
                [row["selected_count"] for row in subset], 50
            ),
        }
    return rows, {
        "seed": SEED,
        "purpose": "computational scaling only; no optimality claim",
        "parameter_count": parameter_count,
        "observed_count": observed_count,
        "repeats_per_pool": repeats,
        "by_pool_size": grouped,
    }


def write_csv(path: Path, rows: list[dict], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    rng = np.random.default_rng(SEED)
    random_rows, random_summary = run_randomized_suite(rng)
    scaling_rows, scaling_summary = run_runtime_scaling(rng)

    write_csv(
        ROOT / "01b_randomized_support_benchmark.csv",
        random_rows,
        [
            "instance_id",
            "dimensions",
            "parameter_count",
            "observed_count",
            "candidate_count",
            "requested_contrasts",
            "cost_sigma",
            "exact_feasible",
            "greedy_feasible",
            "exact_cost",
            "greedy_cost",
            "greedy_cost_gap_fraction",
            "exact_runtime_s",
            "greedy_runtime_s",
            "exact_subsets_evaluated",
            "greedy_candidates_evaluated",
            "exact_selected_count",
            "greedy_selected_count",
        ],
    )
    write_csv(
        ROOT / "01b_greedy_runtime_scaling.csv",
        scaling_rows,
        [
            "pool_size",
            "repeat",
            "feasible",
            "selected_count",
            "candidates_evaluated",
            "runtime_s",
        ],
    )
    (ROOT / "01b_randomized_support_benchmark_full.json").write_text(
        json.dumps(
            {
                "protocol": "V12-20260730", "frozen_parent_protocol": "V10-20260730",
                "randomized_instances": random_rows,
                "runtime_scaling": scaling_rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    summary = {
        "protocol": "V12-20260730", "frozen_parent_protocol": "V10-20260730",
        "randomized_finite_pool": random_summary,
        "greedy_runtime_scaling": scaling_summary,
    }
    (ROOT / "01b_randomized_support_benchmark_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
