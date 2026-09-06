# Auxiliary robust-union analysis, not a source of the main CFD figures.
"""Derive epsilon-coupling costs from the graph-constrained CFD profiles."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUT = HERE / "28_ROBUST_UNION_CFD_PROFILES.json"
CSV_OUTPUT = HERE / "29_ROBUST_UNION_EPSILON_COSTS.csv"
JSON_OUTPUT = HERE / "29_ROBUST_UNION_EPSILON_COSTS.json"
EPSILONS = (0.0, 1.0e-12, 1.0e-10, 1.0e-9, 1.0e-8, 1.0e-7)


def first_attaining(profile, metric, target):
    return min(
        row["relaxation_budget"]
        for row in profile
        if row[metric] >= target
    )


def main():
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = []
    for result in payload["results"]:
        profile = result["profile"]
        full = profile[-1]["euclidean_finite_set_diameter"]
        at_two = next(
            row["euclidean_finite_set_diameter"]
            for row in profile
            if row["relaxation_budget"] == 2.0
        )
        row = {
            "case_id": result["case_id"],
            "support": result["support"],
            "relative_diameter_gap_at_budget_2": (full - at_two) / full,
        }
        for epsilon in EPSILONS:
            key = f"minimum_cost_epsilon_{epsilon:.0e}"
            # The 1e-15 absolute slack only prevents a stored number from
            # failing comparison with itself after JSON decimal round-tripping.
            row[key] = first_attaining(
                profile,
                "euclidean_finite_set_diameter",
                (1.0 - epsilon) * full - 1.0e-15,
            )
        rows.append(row)

    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    output = {
        "definition": (
            "rho_D(epsilon) is the minimum graph-relaxation budget c for which "
            "D(U_c) >= (1-epsilon) D(U_ind)."
        ),
        "epsilons": EPSILONS,
        "rows": rows,
        "interpretation": (
            "Exact attainment and practically certified near-attainment are "
            "reported separately; no epsilon is selected post hoc as universal."
        ),
    }
    JSON_OUTPUT.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
