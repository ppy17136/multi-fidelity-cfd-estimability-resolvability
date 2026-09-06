"""Verification benchmark for finite graph-constrained coupling MILP instances.

Families:
1. tree-structured tied-extrema problems (forest DP versus MILP);
2. general random tied-extrema graphs;
3. three-terminal Multiway-Cut reductions;
4. threshold-crossing scalar decision problems.

These instances check correctness and finite-instance feasibility, not general-graph scalability.
"""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path
import time

import numpy as np

from graph_constrained_coupling import (
    pathway_pair_compatibility,
    tree_graph_relaxation,
)
from graph_coupling_milp import (
    milp_decision_destroying_relaxation,
    milp_directional_full_width,
)


HERE = Path(__file__).resolve().parent
CSV_OUTPUT = HERE / "30_GRAPH_COUPLING_MILP_BENCHMARK.csv"
JSON_OUTPUT = HERE / "30_GRAPH_COUPLING_MILP_BENCHMARK.json"
TIME_LIMIT = 10.0


def elapsed(call):
    start = time.perf_counter()
    value = call()
    return value, time.perf_counter() - start


def random_tree_edges(q, rng):
    return [
        (int(rng.integers(0, child)), child, float(rng.integers(1, 11)))
        for child in range(1, q)
    ]


def random_graph_edges(q, density, rng):
    # Include a random spanning tree so the zero-cost decision baseline is
    # unambiguously one coherent pathway over the whole graph.
    tree = random_tree_edges(q, rng)
    tree_pairs = {tuple(sorted((left, right))) for left, right, _ in tree}
    extra = [
        (left, right, float(rng.integers(1, 11)))
        for left, right in itertools.combinations(range(q), 2)
        if (left, right) not in tree_pairs and rng.random() < density
    ]
    return tree + extra


def tied_projection(q, rng):
    # Integer values deliberately create non-singleton argmax/argmin sets.
    return rng.integers(0, 4, size=(q, 4)).astype(float)


def multiway_allowed(q):
    if q < 3:
        raise ValueError("multiway construction needs three terminals")
    z = np.zeros((q, 4), dtype=float)
    for terminal in range(3):
        z[terminal, 1:] = 0.5
        z[terminal, terminal + 1] = 1.0
    z[3:, 1:] = 1.0
    return pathway_pair_compatibility(z)


def decision_contributions(q, rng):
    if q % 2:
        raise ValueError("paired construction requires even q")
    half = rng.integers(0, 11, size=q // 2).astype(float)
    first = np.concatenate((half, 10.0 - half))
    second = 10.0 - first
    z = np.column_stack((first, second))
    coherent_lower = float(min(z[:, 0].sum(), z[:, 1].sum()))
    independent_lower = float(np.min(z, axis=1).sum())
    threshold = 0.5 * (coherent_lower + independent_lower)
    return z, threshold


def record(family, q, seed, edge_count, certificate, runtime, **extra):
    return {
        "family": family,
        "cells": q,
        "seed": seed,
        "edges": edge_count,
        "status": certificate.status,
        "objective": certificate.objective,
        "dual_bound": certificate.dual_bound,
        "relative_gap": certificate.relative_gap,
        "branch_nodes": certificate.node_count,
        "runtime_s": runtime,
        **extra,
    }


def main():
    rows = []

    for q in (16, 32, 64, 128):
        for seed in range(5):
            rng = np.random.default_rng(2026090400 + 1000 * q + seed)
            allowed = pathway_pair_compatibility(tied_projection(q, rng))
            edges = random_tree_edges(q, rng)
            tree, tree_time = elapsed(lambda: tree_graph_relaxation(allowed, edges))
            cert, runtime = elapsed(
                lambda: milp_directional_full_width(allowed, edges, time_limit=TIME_LIMIT)
            )
            rows.append(record(
                "tree_tied_extrema", q, seed, len(edges), cert, runtime,
                tree_objective=tree.cost,
                tree_runtime_s=tree_time,
                objective_check_abs=(None if cert.objective is None else abs(cert.objective - tree.cost)),
                threshold=None,
                attained_value=None,
            ))

    for q in (12, 20, 32, 48):
        for seed in range(5):
            rng = np.random.default_rng(2026090500 + 1000 * q + seed)
            allowed = pathway_pair_compatibility(tied_projection(q, rng))
            edges = random_graph_edges(q, 0.2, rng)
            cert, runtime = elapsed(
                lambda: milp_directional_full_width(allowed, edges, time_limit=TIME_LIMIT)
            )
            rows.append(record(
                "general_tied_extrema", q, seed, len(edges), cert, runtime,
                tree_objective=None,
                tree_runtime_s=None,
                objective_check_abs=None,
                threshold=None,
                attained_value=None,
            ))

    for q in (15, 25, 40, 60):
        for seed in range(5):
            rng = np.random.default_rng(2026090600 + 1000 * q + seed)
            allowed = multiway_allowed(q)
            edges = random_graph_edges(q, 0.18, rng)
            cert, runtime = elapsed(
                lambda: milp_directional_full_width(allowed, edges, time_limit=TIME_LIMIT)
            )
            rows.append(record(
                "three_terminal_multiway_cut", q, seed, len(edges), cert, runtime,
                tree_objective=None,
                tree_runtime_s=None,
                objective_check_abs=None,
                threshold=None,
                attained_value=None,
            ))

    for q in (16, 24, 32, 48):
        for seed in range(5):
            rng = np.random.default_rng(2026090700 + 1000 * q + seed)
            z, threshold = decision_contributions(q, rng)
            edges = random_graph_edges(q, 0.2, rng)
            (state, cert), runtime = elapsed(
                lambda: milp_decision_destroying_relaxation(
                    z, edges, threshold, time_limit=TIME_LIMIT
                )
            )
            rows.append(record(
                "decision_threshold_crossing", q, seed, len(edges), cert, runtime,
                tree_objective=None,
                tree_runtime_s=None,
                objective_check_abs=None,
                coherent_state=state,
                threshold=threshold,
                attained_value=cert.attained_value,
            ))

    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=sorted(set().union(*(row.keys() for row in rows))))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "graph_coupling_milp_verification_benchmark",
        "time_limit_s": TIME_LIMIT,
        "case_count": len(rows),
        "optimal_count": sum(row["status"] == "optimal" for row in rows),
        "timeout_or_other_count": sum(row["status"] != "optimal" for row in rows),
        "maximum_tree_objective_check_abs": max(
            row["objective_check_abs"] or 0.0
            for row in rows
            if row["family"] == "tree_tied_extrema"
        ),
        "rows": rows,
        "scope": (
            "These instances support correctness and finite-instance feasibility checks only; they do not establish worst-case or general-graph scalability."
        ),
    }
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({key: payload[key] for key in (
        "case_count", "optimal_count", "timeout_or_other_count",
        "maximum_tree_objective_check_abs"
    )}, indent=2))
    for family in sorted(set(row["family"] for row in rows)):
        subset = [row for row in rows if row["family"] == family]
        print(family, {
            "count": len(subset),
            "optimal": sum(row["status"] == "optimal" for row in subset),
            "max_runtime_s": max(row["runtime_s"] for row in subset),
            "max_nodes": max((row["branch_nodes"] or 0) for row in subset),
        })


if __name__ == "__main__":
    main()
