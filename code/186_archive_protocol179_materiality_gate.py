#!/usr/bin/env python3
"""Archive the pre-fit Protocol-179 materiality gate and stop if it fails."""

from __future__ import annotations

import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "186_run_protocol179_fixed_model_comparison.py"
PROTOCOL = ROOT / "184_frozen_protocol179_model_comparison.json"
IMPLEMENTATION = ROOT / "184a_frozen_interaction_evaluation_implementation.json"
OUT = ROOT / "186_protocol179_prefit_materiality_audit.json"
OUT_MD = ROOT / "186_protocol179_prefit_materiality_audit.md"
MARKER_FAIL = ROOT / "PROTOCOL_179_INTERACTION_BELOW_MATERIALITY_FLOOR.ok"
MARKER_REFIT = ROOT / "PROTOCOL_179_MATERIALITY_PERMITS_REFIT.ok"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    for marker in (MARKER_FAIL, MARKER_REFIT):
        if marker.exists():
            marker.unlink()
    spec = importlib.util.spec_from_file_location("runner186", RUNNER)
    if spec is None or spec.loader is None:
        raise RuntimeError("cannot load frozen Protocol-179 runner")
    runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(runner)
    targets = runner.build_exact_targets()

    exact_energy = 0.0
    floor_energy = 0.0
    total_volume = 0.0
    rows = []
    for case_id in runner.PRIMARY_GEOMETRIES:
        target = targets[case_id]
        weight = target["weight"]
        volume = float(np.sum(weight))
        exact_rms = runner.weighted_rms(target["exact_velocity"], weight)
        floor_rms = float(target["mapping_floor_velocity_rms"])
        ratio = exact_rms / floor_rms
        exact_energy += volume * exact_rms**2
        floor_energy += volume * floor_rms**2
        total_volume += volume
        rows.append(
            {
                "case_id": case_id,
                "support_points": len(weight),
                "finite_support_fraction": target[
                    "finite_support_fraction"
                ],
                "exact_velocity_interaction_rms": exact_rms,
                "mapping_floor_velocity_rms": floor_rms,
                "exact_to_mapping_floor_ratio": ratio,
                "three_times_floor_passed": ratio >= 3.0,
                "exact_reattachment_contrast_over_H": target[
                    "exact_reattachment_contrast"
                ],
                "target_sha256": target["target_sha256"],
            }
        )
    pooled_exact = math.sqrt(exact_energy / total_volume)
    pooled_floor = math.sqrt(floor_energy / total_volume)
    pooled_ratio = pooled_exact / pooled_floor
    passed = pooled_ratio >= 3.0
    decision = (
        "materiality_permits_one_fixed_A_C_refit"
        if passed
        else "interaction_not_resolvable_above_numerical_floor_no_refit"
    )
    payload = {
        "protocol": PROTOCOL.name,
        "protocol_sha256": sha256(PROTOCOL),
        "implementation": IMPLEMENTATION.name,
        "implementation_sha256": sha256(IMPLEMENTATION),
        "stage": "before_model_refit_and_before_validation_response_access",
        "per_geometry": rows,
        "pooled": {
            "exact_velocity_interaction_rms": pooled_exact,
            "mapping_floor_velocity_rms": pooled_floor,
            "exact_to_mapping_floor_ratio": pooled_ratio,
            "required_ratio": 3.0,
            "materiality_gate_passed": passed,
        },
        "model_refit_executed": False,
        "validation_response_read": False,
        "test_response_read": False,
        "decision": decision,
        "interpretation": (
            "The exact mesh-closure interaction contrast is too close to the "
            "predeclared numerical mapping floor for a physical attribution "
            "claim. The allowed refit is therefore not executed, because no "
            "model outcome can rescue the failed compound primary endpoint."
        ),
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Protocol 179 pre-fit materiality decision\n\n"
        f"- Pooled exact interaction RMS: `{pooled_exact:.10g}`.\n"
        f"- Pooled frozen mapping floor RMS: `{pooled_floor:.10g}`.\n"
        f"- Signal/floor ratio: `{pooled_ratio:.6f}`; required: `>= 3`.\n"
        f"- Materiality gate: **{'PASS' if passed else 'FAIL'}**.\n"
        f"- Decision: **{decision}**.\n"
        "- No model refit was executed, no validation response was read, and "
        "no sealed test response was opened.\n",
        encoding="utf-8",
    )
    marker = MARKER_REFIT if passed else MARKER_FAIL
    marker.write_text(
        f"decision={decision}\n"
        f"pooled_exact_to_mapping_floor_ratio={pooled_ratio:.12g}\n"
        "validation response read: no\n"
        "test response read: no\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()
