"""Build the compact manuscript table from complete common-partition profiles."""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUT = HERE / "37_REAL_CFD_COMMON_PARTITION_PROFILES.json"
CSV_OUTPUT = HERE / "42_REAL_CFD_COMMON_PARTITION_SUMMARY.csv"
JSON_OUTPUT = HERE / "42_REAL_CFD_COMMON_PARTITION_SUMMARY.json"
BALANCED_BUDGET = 0.5
EPSILON_KEY = "1e-07"


def at_budget(profile, budget):
    matches = [row for row in profile if abs(row["relaxation_budget"] - budget) < 1e-12]
    if len(matches) != 1:
        raise RuntimeError(f"expected one row at budget {budget}, found {len(matches)}")
    return matches[0]


def main():
    source = json.loads(INPUT.read_text(encoding="utf-8"))
    rows = []
    for audit in source["results"]:
        balanced = next(
            model for model in audit["cost_models"]
            if model["cost_model"] == "balanced"
        )
        zero = at_budget(balanced["profile"], 0.0)
        half = at_budget(balanced["profile"], BALANCED_BUDGET)
        independent = audit["independent_metrics"]
        rows.append({
            "case_id": audit["case_id"],
            "support": audit["support"],
            "R_zero_cost": zero["minimum_single_partition_ratio"],
            "R_balanced_cost_0p5": half["minimum_single_partition_ratio"],
            "R_independent_endpoint": independent["ratio"],
            "D_zero_cost": zero["maximum_single_partition_diameter"],
            "D_balanced_cost_0p5": half["maximum_single_partition_diameter"],
            "D_independent_endpoint": independent["diameter"],
            "relative_D_gap_at_0p5": half["diameter_relative_gap"],
            "exact_independent_D_cost": balanced["minimum_cost_exact_independent_diameter"],
            "D_cost_at_relative_tolerance_1e-7": balanced["minimum_diameter_cost_by_epsilon"][EPSILON_KEY],
            "D_witness_partition_at_0p5": half["diameter_partition"],
            "R_witness_partition_at_0p5": half["ratio_partition"],
        })

    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    zero_min = min(row["R_zero_cost"] for row in rows)
    half_max = max(row["R_balanced_cost_0p5"] for row in rows)
    payload = {
        "status": "common_partition_manuscript_summary",
        "source": INPUT.name,
        "audit_count": len(rows),
        "balanced_model": {
            "edge_cost": 0.25,
            "total_graph_cost": 1.0,
            "reported_budget": BALANCED_BUDGET,
        },
        "all_audit_common_threshold_interval": {
            "lower_open": half_max,
            "upper_closed": zero_min,
            "nonempty": half_max < zero_min,
        },
        "rows": rows,
        "claim_boundary": (
            "Each S, D, and R value is computed within one common component "
            "partition. Edge costs are declared normalized sensitivity "
            "weights, not fitted probabilities or physical error scales."
        ),
    }
    JSON_OUTPUT.write_bytes(json.dumps(payload, indent=2).encode("utf-8"))
    print(json.dumps({
        "audit_count": len(rows),
        "threshold_interval": payload["all_audit_common_threshold_interval"],
        "exact_cost_counts": {
            str(cost): sum(row["exact_independent_D_cost"] == cost for row in rows)
            for cost in sorted({row["exact_independent_D_cost"] for row in rows})
        },
        "epsilon_costs": sorted({
            row["D_cost_at_relative_tolerance_1e-7"] for row in rows
        }),
    }, indent=2))


if __name__ == "__main__":
    main()
