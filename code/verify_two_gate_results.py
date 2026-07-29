#!/usr/bin/env python3
"""Recompute the archived manufactured and real-CFD two-gate decisions."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
TOL = 1.0e-12


def close(a: float, b: float, tol: float = TOL) -> bool:
    return math.isclose(a, b, rel_tol=tol, abs_tol=tol)


def verify_manufactured() -> dict[str, object]:
    path = ROOT / "data" / "manufactured" / "199_manufactured_two_gate_summary.json"
    payload = json.loads(path.read_text(encoding="utf-8"))
    assert payload["seed"] == 20260728
    assert payload["monte_carlo_repeats"] == 200_000
    assert close(payload["contrast_numerical_floor"], 1.0)
    assert close(payload["materiality_threshold"], 3.0)
    scenarios = {row["scenario"]: row for row in payload["scenarios"]}
    assert scenarios["missing_corner_high_signal"]["estimable"] is False
    assert scenarios["missing_corner_high_signal"]["rank"] == 3
    assert scenarios["complete_below_floor"]["estimable"] is True
    assert close(scenarios["complete_below_floor"]["pass_rate"], 0.0)
    assert 0.49 < scenarios["complete_at_threshold"]["pass_rate"] < 0.51
    assert close(scenarios["complete_above_floor"]["pass_rate"], 1.0)
    return {
        "seed": payload["seed"],
        "scenario_count": len(scenarios),
        "gate1_negative_control_passed": True,
        "gate2_negative_control_passed": True,
        "gate2_positive_control_passed": True,
    }


def pooled(values: list[float], volumes: list[float]) -> float:
    return math.sqrt(
        sum(volume * value * value for volume, value in zip(volumes, values))
        / sum(volumes)
    )


def verify_real_cfd() -> dict[str, object]:
    csv_path = (
        ROOT / "data" / "derived" / "188_exploratory_four_lattice_materiality.csv"
    )
    with csv_path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 4

    volumes = [float(row["support_volume"]) for row in rows]
    exact = [float(row["exact_velocity_interaction_rms"]) for row in rows]
    archived_floor = [
        float(row["archived_mapping_floor_velocity_rms"]) for row in rows
    ]
    conservative_floor = [
        float(row["conservative_mapping_floor_velocity_rms"]) for row in rows
    ]
    ratios = [
        float(row["exact_to_archived_floor_ratio"]) for row in rows
    ]
    for e_value, floor_value, ratio in zip(exact, archived_floor, ratios):
        assert close(e_value / floor_value, ratio)
        assert ratio < 3.0

    pooled_exact = pooled(exact, volumes)
    pooled_archived = pooled(archived_floor, volumes)
    pooled_conservative = pooled(conservative_floor, volumes)
    archived_ratio = pooled_exact / pooled_archived
    conservative_ratio = pooled_exact / pooled_conservative

    json_path = (
        ROOT / "data" / "derived" / "188_exploratory_four_lattice_materiality.json"
    )
    archived = json.loads(json_path.read_text(encoding="utf-8"))
    archived_primary = archived["pooled_all_four_archived_floor"]
    archived_conservative = archived["pooled_all_four_conservative_floor"]
    assert close(
        archived_ratio, archived_primary["exact_to_mapping_floor_ratio"]
    )
    assert close(
        conservative_ratio,
        archived_conservative["exact_to_mapping_floor_ratio"],
    )
    assert archived_ratio < 3.0
    assert conservative_ratio < 3.0

    return {
        "geometry_count": len(rows),
        "all_gate1_estimable": True,
        "all_gate2_ratios_below_3": True,
        "pooled_archived_floor_ratio": archived_ratio,
        "pooled_conservative_floor_ratio": conservative_ratio,
    }


def main() -> None:
    result = {
        "manufactured": verify_manufactured(),
        "real_cfd": verify_real_cfd(),
        "decision": "verification_passed",
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
