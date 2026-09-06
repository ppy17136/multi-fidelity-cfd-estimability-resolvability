# Auxiliary robust-union analysis, not a source of the main CFD figures.
"""Compute graph-constrained coupling profiles for eight released CFD audits.

The four fidelity cells form the declared two-axis square

    F000 --mesh-- F100
      |              |
    closure        closure
      |              |
    F010 --mesh-- F110

Every edge has unit cost in this first, deliberately non-empirical reference
model. A pathway tuple pays one unit on an edge when its endpoints use
different mapping pathways. Results are written beside this script as CSV
and JSON.
"""

from __future__ import annotations

import csv
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE.parents[1] / "data" / "cfd_mapping"
CSV_OUTPUT = HERE / "28_ROBUST_UNION_CFD_PROFILES.csv"
JSON_OUTPUT = HERE / "28_ROBUST_UNION_CFD_PROFILES.json"
RTOL = 1.0e-10


def square_edges(fidelities):
    index = {name: i for i, name in enumerate(fidelities)}
    expected = {"F000", "F100", "F010", "F110"}
    if set(index) != expected:
        raise ValueError(f"expected the four two-axis fidelity cells, got {set(index)}")
    return (
        (index["F000"], index["F100"], "mesh", 1.0),
        (index["F010"], index["F110"], "mesh", 1.0),
        (index["F000"], index["F010"], "closure", 1.0),
        (index["F100"], index["F110"], "closure", 1.0),
    )


def exact_metrics(vectors):
    norms = np.linalg.norm(vectors, axis=1)
    gram = vectors @ vectors.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    diameter = float(np.sqrt(float(distance2.max())))
    signal = float(norms.min())
    ratio = float(signal / diameter) if diameter > 0 else float("inf")
    return signal, diameter, ratio


def evaluate(path):
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(x) for x in data["pathways"]]
        fidelities = [str(x) for x in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])

    pathway_count = len(pathways)
    q = len(fidelities)
    contributions = np.stack([
        np.stack([
            coefficients[i] * raw[pathway * q + i]
            for pathway in range(pathway_count)
        ])
        for i in range(q)
    ])
    choices = list(itertools.product(range(pathway_count), repeat=q))
    vectors = np.stack([
        sum(contributions[i, choice[i]] for i in range(q))
        for choice in choices
    ])
    edges = square_edges(fidelities)
    costs = np.asarray([
        sum(weight for left, right, _, weight in edges if choice[left] != choice[right])
        for choice in choices
    ])
    budgets = sorted(set(float(value) for value in costs))
    full_signal, full_diameter, _ = exact_metrics(vectors)
    rows = []
    for budget in budgets:
        selected = costs <= budget + 1.0e-12
        signal, diameter, ratio = exact_metrics(vectors[selected])
        rows.append({
            "case_id": case_id,
            "support": support,
            "relaxation_budget": budget,
            "admissible_members": int(np.count_nonzero(selected)),
            "minimum_signal": signal,
            "euclidean_finite_set_diameter": diameter,
            "signal_to_diameter_ratio": ratio,
            "attains_independent_minimum_signal": bool(np.isclose(signal, full_signal, rtol=RTOL, atol=1.0e-14)),
            "attains_independent_diameter": bool(np.isclose(diameter, full_diameter, rtol=RTOL, atol=1.0e-14)),
        })
    diameter_cost = min(row["relaxation_budget"] for row in rows if row["attains_independent_diameter"])
    signal_cost = min(row["relaxation_budget"] for row in rows if row["attains_independent_minimum_signal"])
    return {
        "case_id": case_id,
        "support": support,
        "fidelities": fidelities,
        "pathways": pathways,
        "edge_model": [
            {"left": fidelities[left], "right": fidelities[right], "axis": axis, "cost": weight}
            for left, right, axis, weight in edges
        ],
        "minimum_unit_edge_cut_cost_for_independent_diameter": diameter_cost,
        "minimum_unit_edge_cut_cost_for_independent_minimum_signal": signal_cost,
        "profile": rows,
    }


def main():
    results = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    rows = [row for result in results for row in result["profile"]]
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "graph_constrained_coupling_profile",
        "graph": "two-axis fidelity square",
        "edge_cost_model": "unit cost per released mesh- or closure-axis pathway disagreement",
        "cost_interpretation": (
            "Dimensionless reference costs expose topology; they are not "
            "empirically identified probabilities or physical error magnitudes."
        ),
        "audit_count": len(results),
        "all_diameters_attained_by_budget_two_or_less": all(
            result["minimum_unit_edge_cut_cost_for_independent_diameter"] <= 2.0
            for result in results
        ),
        "results": results,
    }
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "csv": str(CSV_OUTPUT),
        "json": str(JSON_OUTPUT),
        "audit_count": len(results),
        "diameter_costs": [result["minimum_unit_edge_cut_cost_for_independent_diameter"] for result in results],
        "signal_costs": [result["minimum_unit_edge_cut_cost_for_independent_minimum_signal"] for result in results],
    }, indent=2))


if __name__ == "__main__":
    main()
