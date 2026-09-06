"""Exact common-partition graph profiles for the eight CFD mapping audits.

This script distinguishes two valid but different constructions:

1. a robust union of members attainable under any relaxation up to a budget;
2. one common removed-edge set (partition) whose uncertainty set must support
   all members used in a set functional.

The graph-labelling theorem concerns construction 2. Every edge-removal mask
of the four-cell square is enumerated, duplicate component partitions retain
their minimum cut cost, and metrics are computed within each single partition.
"""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "data" / "cfd_mapping"
CSV_OUTPUT = HERE / "37_REAL_CFD_COMMON_PARTITION_PROFILES.csv"
JSON_OUTPUT = HERE / "37_REAL_CFD_COMMON_PARTITION_PROFILES.json"
EPSILONS = (0.0, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7)
COST_MODELS = {
    "balanced": {"mesh": 0.25, "closure": 0.25},
    "mesh_dominant": {"mesh": 0.40, "closure": 0.10},
    "closure_dominant": {"mesh": 0.10, "closure": 0.40},
}


def canonical_partition(q, retained_edges):
    parent = list(range(q))

    def find(node):
        while parent[node] != node:
            parent[node] = parent[parent[node]]
            node = parent[node]
        return node

    def union(left, right):
        left_root, right_root = find(left), find(right)
        if left_root != right_root:
            parent[right_root] = left_root

    for left, right in retained_edges:
        union(left, right)
    groups = {}
    for node in range(q):
        groups.setdefault(find(node), []).append(node)
    return tuple(sorted(tuple(group) for group in groups.values()))


def partitions_with_minimum_cost(q, edges, axis_cost):
    best = {}
    for mask in range(1 << len(edges)):
        removed = [bool(mask & (1 << edge)) for edge in range(len(edges))]
        retained = [(left, right) for remove, (left, right, _) in zip(removed, edges)
                    if not remove]
        partition = canonical_partition(q, retained)
        cost = sum(axis_cost[axis] for remove, (_, _, axis) in zip(removed, edges)
                   if remove)
        if partition not in best or cost < best[partition]:
            best[partition] = cost
    return tuple(sorted(best.items(), key=lambda item: (item[1], item[0])))


def allowed_indices(partition, pathway_count, choice_to_index):
    indices = []
    for component_labels in itertools.product(range(pathway_count), repeat=len(partition)):
        choice = [None] * sum(map(len, partition))
        for component, pathway in zip(partition, component_labels):
            for cell in component:
                choice[cell] = pathway
        indices.append(choice_to_index[tuple(choice)])
    return np.asarray(indices, dtype=int)


def metrics(norms, distance2, indices):
    signal = float(norms[indices].min())
    diameter = float(np.sqrt(distance2[np.ix_(indices, indices)].max()))
    return signal, diameter, signal / diameter if diameter > 0.0 else float("inf")


def evaluate(path):
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(value) for value in data["pathways"]]
        fidelities = [str(value) for value in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])
    p, q = len(pathways), len(fidelities)
    contributions = np.stack([
        np.stack([coefficients[cell] * raw[pathway * q + cell]
                  for pathway in range(p)])
        for cell in range(q)
    ])
    choices = list(itertools.product(range(p), repeat=q))
    choice_to_index = {choice: index for index, choice in enumerate(choices)}
    vectors = np.stack([
        sum(contributions[cell, choice[cell]] for cell in range(q))
        for choice in choices
    ])
    norms = np.linalg.norm(vectors, axis=1)
    gram = vectors @ vectors.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    full_indices = np.arange(len(choices))
    full_signal, full_diameter, full_ratio = metrics(norms, distance2, full_indices)
    cell = {name: index for index, name in enumerate(fidelities)}
    edges = (
        (cell["F000"], cell["F100"], "mesh"),
        (cell["F010"], cell["F110"], "mesh"),
        (cell["F000"], cell["F010"], "closure"),
        (cell["F100"], cell["F110"], "closure"),
    )

    model_results = []
    rows = []
    for model_name, axis_cost in COST_MODELS.items():
        records = []
        for partition, cost in partitions_with_minimum_cost(q, edges, axis_cost):
            indices = allowed_indices(partition, p, choice_to_index)
            signal, diameter, ratio = metrics(norms, distance2, indices)
            records.append({
                "partition": partition,
                "cost": float(cost),
                "members": int(len(indices)),
                "signal": signal,
                "diameter": diameter,
                "ratio": ratio,
            })
        budgets = sorted(set(record["cost"] for record in records))
        profile = []
        for budget in budgets:
            admissible = [record for record in records if record["cost"] <= budget + 1.0e-12]
            best_diameter = max(admissible, key=lambda record: record["diameter"])
            worst_ratio = min(admissible, key=lambda record: record["ratio"])
            minimum_signal = min(admissible, key=lambda record: record["signal"])
            row = {
                "case_id": case_id,
                "support": support,
                "cost_model": model_name,
                "relaxation_budget": budget,
                "admissible_component_partitions": len(admissible),
                "maximum_single_partition_diameter": best_diameter["diameter"],
                "diameter_partition": json.dumps(best_diameter["partition"]),
                "minimum_single_partition_signal": minimum_signal["signal"],
                "minimum_single_partition_ratio": worst_ratio["ratio"],
                "ratio_partition": json.dumps(worst_ratio["partition"]),
                "diameter_relative_gap": (full_diameter - best_diameter["diameter"]) / full_diameter,
            }
            rows.append(row)
            profile.append(row)
        exact_cost = min(
            row["relaxation_budget"] for row in profile
            if row["maximum_single_partition_diameter"] >= full_diameter
        )
        epsilon_costs = {
            f"{epsilon:.0e}": min(
                row["relaxation_budget"] for row in profile
                if row["maximum_single_partition_diameter"] >= (1.0 - epsilon) * full_diameter
            )
            for epsilon in EPSILONS
        }
        model_results.append({
            "cost_model": model_name,
            "axis_edge_cost": axis_cost,
            "distinct_component_partitions": len(records),
            "minimum_cost_exact_independent_diameter": exact_cost,
            "minimum_diameter_cost_by_epsilon": epsilon_costs,
            "profile": profile,
        })
    return {
        "case_id": case_id,
        "support": support,
        "independent_metrics": {
            "signal": full_signal,
            "diameter": full_diameter,
            "ratio": full_ratio,
        },
        "cost_models": model_results,
    }, rows


def main():
    packed = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    results = [item[0] for item in packed]
    rows = [row for item in packed for row in item[1]]
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "exact_single_common_partition_profiles",
        "cost_models": COST_MODELS,
        "audit_count": len(results),
        "results": results,
        "claim_boundary": (
            "Every reported set functional is attained within one common "
            "component partition. Edge weights are declared sensitivity models, "
            "not empirically identified probabilities or physical errors."
        ),
    }
    JSON_OUTPUT.write_bytes(json.dumps(payload, indent=2).encode("utf-8"))
    compact = {}
    for model_name in COST_MODELS:
        compact[model_name] = {
            "exact_diameter_costs": [
                next(model for model in result["cost_models"]
                     if model["cost_model"] == model_name)
                ["minimum_cost_exact_independent_diameter"]
                for result in results
            ],
            "epsilon_1e-7_costs": [
                next(model for model in result["cost_models"]
                     if model["cost_model"] == model_name)
                ["minimum_diameter_cost_by_epsilon"]["1e-07"]
                for result in results
            ],
        }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
