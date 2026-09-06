# Auxiliary robust-union analysis, not a source of the main CFD figures.
"""Build threshold-to-minimum-relaxation frontiers for the eight CFD audits.

For a relaxation budget c, U_c is nested in c. Hence its minimum norm S(c)
is nonincreasing, its finite-set diameter D(c) is nondecreasing, and
R(c)=S(c)/D(c) is nonincreasing. If the zero-cost state passes a threshold
tau under the rule R >= tau, the first budget with R(c) < tau is an exact
minimum robust-union ratio-classification-change budget for the declared graph.
This is not the scalar decision certificate or a common-partition profile.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path


HERE = Path(__file__).resolve().parent
INPUT = HERE / "28_ROBUST_UNION_CFD_PROFILES.csv"
CSV_OUTPUT = HERE / "33_ROBUST_UNION_THRESHOLD_FRONTIERS.csv"
JSON_OUTPUT = HERE / "33_ROBUST_UNION_THRESHOLD_FRONTIERS.json"


def threshold_intervals(profile):
    """Intervals of tau having each positive minimum destroying budget.

    With pass defined by R >= tau, budget c first destroys a zero-cost pass
    for R(c) < tau <= min(previous attained ratios). Open lower and closed
    upper endpoints retain the strict failure convention.
    """
    zero_ratio = profile[0]["ratio"]
    incumbent = zero_ratio
    intervals = []
    for row in profile[1:]:
        ratio = row["ratio"]
        if ratio < incumbent:
            intervals.append({
                "minimum_destroying_budget": row["budget"],
                "threshold_lower_open": ratio,
                "threshold_upper_closed": incumbent,
            })
            incumbent = ratio
    return intervals


def main():
    grouped = {}
    with INPUT.open(newline="", encoding="utf-8") as handle:
        for raw in csv.DictReader(handle):
            key = (raw["case_id"], raw["support"])
            grouped.setdefault(key, []).append({
                "budget": float(raw["relaxation_budget"]),
                "ratio": float(raw["signal_to_diameter_ratio"]),
                "signal": float(raw["minimum_signal"]),
                "diameter": float(raw["euclidean_finite_set_diameter"]),
            })

    results = []
    flat = []
    for (case_id, support), profile in sorted(grouped.items()):
        profile.sort(key=lambda row: row["budget"])
        if any(profile[i + 1]["ratio"] > profile[i]["ratio"] for i in range(len(profile) - 1)):
            raise AssertionError("robustness ratio must be nonincreasing")
        intervals = threshold_intervals(profile)
        result = {
            "case_id": case_id,
            "support": support,
            "profile": profile,
            "threshold_cost_intervals": intervals,
        }
        results.append(result)
        for interval in intervals:
            flat.append({"case_id": case_id, "support": support, **interval})

    cost_two = [
        next(interval for interval in result["threshold_cost_intervals"]
             if interval["minimum_destroying_budget"] == 2.0)
        for result in results
    ]
    common_lower = max(interval["threshold_lower_open"] for interval in cost_two)
    common_upper = min(interval["threshold_upper_closed"] for interval in cost_two)
    payload = {
        "status": "exact_threshold_to_relaxation_cost_frontier",
        "decision_rule": "pass iff R(c) >= tau; destroy when R(c) < tau",
        "monotonicity": (
            "Nested U_c implies nonincreasing S(c), nondecreasing D(c), "
            "and nonincreasing R(c)=S(c)/D(c)."
        ),
        "audit_count": len(results),
        "common_minimum_cost_two_interval": {
            "lower_open": common_lower,
            "upper_closed": common_upper,
            "nonempty": common_lower < common_upper,
        },
        "results": results,
        "claim_boundary": (
            "The frontier is exact for the declared unit-cost square graph and "
            "released finite mapping sets, not for total CFD uncertainty."
        ),
    }
    with CSV_OUTPUT.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(flat[0]))
        writer.writeheader()
        writer.writerows(flat)
    JSON_OUTPUT.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    print(json.dumps({
        "audit_count": len(results),
        "interval_row_count": len(flat),
        "common_minimum_cost_two_interval": payload["common_minimum_cost_two_interval"],
    }, indent=2))


if __name__ == "__main__":
    main()
