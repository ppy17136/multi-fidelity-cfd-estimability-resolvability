"""Report the intrinsic linear dimension of the released mapping sets.

The affine variation of a q-cell, P-pathway recombined set is contained in
the span of the q(P-1) within-cell pathway differences.  The script is
read-only and prints a JSON report to stdout.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "data" / "cfd_mapping"


def evaluate(path: Path) -> dict:
    with np.load(path) as data:
        raw = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(x) for x in data["pathways"]]
        fidelities = [str(x) for x in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])

    p_count = len(pathways)
    q = len(fidelities)
    contributions = np.stack([
        np.stack([
            coefficients[i] * raw[p * q + i]
            for p in range(p_count)
        ])
        for i in range(q)
    ])

    differences = np.stack([
        contributions[i, p] - contributions[i, 0]
        for i in range(q)
        for p in range(1, p_count)
    ])
    singular_values = np.linalg.svd(differences, compute_uv=False)
    largest = float(singular_values[0]) if singular_values.size else 0.0

    ranks = {}
    for relative_tolerance in (1.0e-8, 1.0e-10, 1.0e-12):
        threshold = relative_tolerance * largest
        ranks[f"rtol_{relative_tolerance:.0e}"] = int(
            np.count_nonzero(singular_values > threshold)
        )

    return {
        "case_id": case_id,
        "support": support,
        "ambient_dimension": int(contributions.shape[-1]),
        "cell_count": q,
        "pathway_count": p_count,
        "difference_generator_count": int(differences.shape[0]),
        "rank_upper_bound_q_times_P_minus_1": q * (p_count - 1),
        "numerical_ranks": ranks,
        "singular_values": [float(x) for x in singular_values],
    }


def main() -> None:
    results = [evaluate(path) for path in sorted(INPUT.glob("*.npz"))]
    payload = {
        "input_count": len(results),
        "all_affine_ranks_at_most_q_times_P_minus_1": all(
            max(item["numerical_ranks"].values())
            <= item["rank_upper_bound_q_times_P_minus_1"]
            for item in results
        ),
        "results": results,
    }
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
