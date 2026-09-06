"""Numerically stable audit of near-tied graph-coupling diameters.

The complete 256-member finite sets are first enumerated with the same
float64 Gram calculation used by the reference profile.  For each budget,
the leading candidate pairs are then recomputed from direct vector
differences and accumulated with :func:`math.fsum`.  The second calculation
avoids cancellation in ``||x||^2 + ||y||^2 - 2 x^T y``.

This is an audit at the precision of the released float64 caches.  It does
not infer accuracy beyond those caches.
"""

from __future__ import annotations

import csv
import itertools
import json
import math
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "data" / "cfd_mapping"
CSV_OUTPUT = HERE / "31_REAL_CFD_DIAMETER_NEAR_TIE_AUDIT.csv"
JSON_OUTPUT = HERE / "31_REAL_CFD_DIAMETER_NEAR_TIE_AUDIT.json"
TOP_CANDIDATES = 12


def stable_distance(left: np.ndarray, right: np.ndarray) -> float:
    difference = left - right
    squared = math.fsum(float(value) * float(value) for value in difference)
    return math.sqrt(squared)


def square_cost(choice, fidelity_index):
    edges = (
        (fidelity_index["F000"], fidelity_index["F100"]),
        (fidelity_index["F010"], fidelity_index["F110"]),
        (fidelity_index["F000"], fidelity_index["F010"]),
        (fidelity_index["F100"], fidelity_index["F110"]),
    )
    return sum(choice[left] != choice[right] for left, right in edges)


def leading_pairs(vectors, member_indices):
    selected = vectors[member_indices]
    gram = selected @ selected.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    upper_left, upper_right = np.triu_indices(len(member_indices), k=1)
    values = distance2[upper_left, upper_right]
    count = min(TOP_CANDIDATES, len(values))
    leading = np.argpartition(values, -count)[-count:]
    pairs = [
        (int(member_indices[upper_left[index]]), int(member_indices[upper_right[index]]))
        for index in leading
    ]
    return pairs, float(np.sqrt(values.max()))


def evaluate(path):
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(value) for value in data["pathways"]]
        fidelities = [str(value) for value in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])

    pathway_count = len(pathways)
    q = len(fidelities)
    contributions = np.stack([
        np.stack([
            coefficients[cell] * raw[pathway * q + cell]
            for pathway in range(pathway_count)
        ])
        for cell in range(q)
    ])
    choices = list(itertools.product(range(pathway_count), repeat=q))
    vectors = np.stack([
        sum(contributions[cell, choice[cell]] for cell in range(q))
        for choice in choices
    ])
    fidelity_index = {name: index for index, name in enumerate(fidelities)}
    costs = np.asarray([square_cost(choice, fidelity_index) for choice in choices])

    audit = {}
    for budget in (2, 4):
        members = np.flatnonzero(costs <= budget)
        candidate_pairs, gram_diameter = leading_pairs(vectors, members)
        stable = [
            (stable_distance(vectors[left], vectors[right]), left, right)
            for left, right in candidate_pairs
        ]
        distance, left, right = max(stable)
        audit[budget] = {
            "diameter_direct_fsum": distance,
            "diameter_gram": gram_diameter,
            "gram_minus_direct": gram_diameter - distance,
            "left_member": left,
            "right_member": right,
            "left_choice": list(choices[left]),
            "right_choice": list(choices[right]),
            "left_cost": int(costs[left]),
            "right_cost": int(costs[right]),
        }

    budget_two = audit[2]["diameter_direct_fsum"]
    full = audit[4]["diameter_direct_fsum"]
    absolute_gap = full - budget_two
    relative_gap = absolute_gap / full
    ulp = math.ulp(full)
    return {
        "case_id": case_id,
        "support": support,
        "vector_dimension": int(vectors.shape[1]),
        "budget_two_diameter_direct_fsum": budget_two,
        "full_diameter_direct_fsum": full,
        "absolute_gap": absolute_gap,
        "relative_gap": relative_gap,
        "full_diameter_ulp": ulp,
        "gap_in_full_diameter_ulps": absolute_gap / ulp,
        "strictly_larger_at_budget_four": bool(absolute_gap > 0.0),
        "budget_two_certificate": audit[2],
        "full_certificate": audit[4],
    }


def main():
    results = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    rows = [{key: value for key, value in result.items() if not key.endswith("_certificate")}
            for result in results]
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    payload = {
        "status": "released_float64_near_tie_audit",
        "method": (
            "Complete Gram enumeration followed by direct-difference math.fsum "
            f"recalculation of the leading {TOP_CANDIDATES} pairs per budget."
        ),
        "claim_boundary": "No accuracy is inferred beyond the released float64 caches.",
        "all_strict_signs_confirmed": all(
            result["strictly_larger_at_budget_four"] for result in results
        ),
        "results": results,
    }
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "audit_count": len(results),
        "all_strict_signs_confirmed": payload["all_strict_signs_confirmed"],
        "relative_gaps": [result["relative_gap"] for result in results],
        "gap_ulps": [result["gap_in_full_diameter_ulps"] for result in results],
    }, indent=2))


if __name__ == "__main__":
    main()
