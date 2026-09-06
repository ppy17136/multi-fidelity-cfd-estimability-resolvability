# Auxiliary robust-union analysis, not a source of the main CFD figures.
"""Sensitivity of CFD graph-coupling profiles to declared axis costs.

All models have total graph weight one. They are illustrative provenance
credibility models, not fitted probabilities or physical-error magnitudes.
"""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE.parents[1] / "data" / "cfd_mapping"
CSV_OUTPUT = HERE / "35_ROBUST_UNION_EDGE_COST_SENSITIVITY.csv"
JSON_OUTPUT = HERE / "35_ROBUST_UNION_EDGE_COST_SENSITIVITY.json"
EPSILONS = (0.0, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7)
COST_MODELS = {
    "balanced": {"mesh": 0.25, "closure": 0.25},
    "mesh_dominant": {"mesh": 0.40, "closure": 0.10},
    "closure_dominant": {"mesh": 0.10, "closure": 0.40},
}


def geometry(vectors):
    norms = np.linalg.norm(vectors, axis=1)
    gram = vectors @ vectors.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    return norms, distance2


def subset_metrics(norms, distance2, selected):
    indices = np.flatnonzero(selected)
    signal = float(norms[indices].min())
    diameter = float(np.sqrt(distance2[np.ix_(indices, indices)].max()))
    ratio = signal / diameter if diameter > 0.0 else float("inf")
    return signal, diameter, ratio


def evaluate(path):
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(value) for value in data["pathways"]]
        fidelities = [str(value) for value in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])
    p = len(pathways)
    q = len(fidelities)
    contributions = np.stack([
        np.stack([coefficients[cell] * raw[pathway * q + cell]
                  for pathway in range(p)])
        for cell in range(q)
    ])
    choices = list(itertools.product(range(p), repeat=q))
    vectors = np.stack([
        sum(contributions[cell, choice[cell]] for cell in range(q))
        for choice in choices
    ])
    norms, distance2 = geometry(vectors)
    full_signal, full_diameter, _ = subset_metrics(
        norms, distance2, np.ones(len(choices), dtype=bool)
    )
    index = {name: cell for cell, name in enumerate(fidelities)}
    edges = (
        (index["F000"], index["F100"], "mesh"),
        (index["F010"], index["F110"], "mesh"),
        (index["F000"], index["F010"], "closure"),
        (index["F100"], index["F110"], "closure"),
    )
    rows = []
    summaries = []
    for model_name, axis_cost in COST_MODELS.items():
        costs = np.asarray([
            sum(axis_cost[axis] for left, right, axis in edges
                if choice[left] != choice[right])
            for choice in choices
        ])
        budgets = sorted(set(round(float(cost), 12) for cost in costs))
        profile = []
        for budget in budgets:
            selected = costs <= budget + 1.0e-12
            signal, diameter, ratio = subset_metrics(norms, distance2, selected)
            row = {
                "case_id": case_id,
                "support": support,
                "cost_model": model_name,
                "mesh_edge_cost": axis_cost["mesh"],
                "closure_edge_cost": axis_cost["closure"],
                "relaxation_budget": budget,
                "admissible_members": int(selected.sum()),
                "minimum_signal": signal,
                "diameter": diameter,
                "ratio": ratio,
                "diameter_relative_gap": (full_diameter - diameter) / full_diameter,
                "signal_relative_gap": (signal - full_signal) / max(abs(full_signal), 1.0e-300),
            }
            rows.append(row)
            profile.append(row)
        exact_cost = min(
            row["relaxation_budget"] for row in profile
            if row["diameter"] >= full_diameter
        )
        signal_cost = min(
            row["relaxation_budget"] for row in profile
            if np.isclose(row["minimum_signal"], full_signal, rtol=1.0e-10, atol=1.0e-14)
        )
        epsilon_costs = {
            f"{epsilon:.0e}": min(
                row["relaxation_budget"] for row in profile
                if row["diameter"] >= (1.0 - epsilon) * full_diameter
            )
            for epsilon in EPSILONS
        }
        summaries.append({
            "cost_model": model_name,
            "minimum_cost_exact_diameter": exact_cost,
            "minimum_cost_independent_signal": signal_cost,
            "minimum_diameter_cost_by_epsilon": epsilon_costs,
        })
    return {"case_id": case_id, "support": support, "summaries": summaries}, rows


def main():
    packed = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    results = [item[0] for item in packed]
    rows = [row for item in packed for row in item[1]]
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "declared_edge_cost_sensitivity",
        "normalization": "Each four-edge graph has total edge weight one.",
        "cost_models": COST_MODELS,
        "epsilon_values": list(EPSILONS),
        "audit_count": len(results),
        "results": results,
        "claim_boundary": (
            "The three cost models are sensitivity probes. Their weights are "
            "not empirically identified probabilities or physical errors."
        ),
    }
    JSON_OUTPUT.write_bytes(json.dumps(payload, indent=2).encode("utf-8"))
    compact = {}
    for model in COST_MODELS:
        compact[model] = {
            "exact_diameter_costs": [
                next(x for x in result["summaries"] if x["cost_model"] == model)
                ["minimum_cost_exact_diameter"]
                for result in results
            ],
            "epsilon_1e-7_costs": [
                next(x for x in result["summaries"] if x["cost_model"] == model)
                ["minimum_diameter_cost_by_epsilon"]["1e-07"]
                for result in results
            ],
        }
    print(json.dumps(compact, indent=2))


if __name__ == "__main__":
    main()
