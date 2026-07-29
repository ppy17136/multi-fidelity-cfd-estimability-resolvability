#!/usr/bin/env python3
"""Finalize the scientific disposition of the Protocol-179 branch."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
INTEGRATION = ROOT / "183_second_stage_integration_audit.json"
FREEZE = ROOT / "184_frozen_protocol179_model_comparison.json"
IMPLEMENTATION = ROOT / "184a_frozen_interaction_evaluation_implementation.json"
INPUTS = ROOT / "185_protocol179_model_input_summary.json"
MATERIALITY = ROOT / "186_protocol179_prefit_materiality_audit.json"
OUT_JSON = ROOT / "187_protocol179_scientific_decision.json"
OUT_MD = ROOT / "187_protocol179_scientific_decision.md"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    integration = json.loads(INTEGRATION.read_text(encoding="utf-8"))
    inputs = json.loads(INPUTS.read_text(encoding="utf-8"))
    materiality = json.loads(MATERIALITY.read_text(encoding="utf-8"))
    pooled = materiality["pooled"]
    decision = {
        "date": "2026-07-28",
        "branch": "Protocol 179 sequential real-CFD confirmation",
        "execution_result": {
            "second_stage_labels_completed": 3,
            "second_stage_labels_platform_accepted": 3,
            "combined_new_platform_accepted_labels": 8,
            "design_A_CFD_labels": 24,
            "design_C_CFD_labels_prepared": 32,
            "strict_complete_lattice_count": integration[
                "strict_complete_lattice_count"
            ],
            "strict_geometry_affine_rank": integration[
                "strict_geometry_affine_rank"
            ],
        },
        "primary_gate_result": {
            "pooled_exact_interaction_rms": pooled[
                "exact_velocity_interaction_rms"
            ],
            "pooled_mapping_floor_rms": pooled[
                "mapping_floor_velocity_rms"
            ],
            "signal_to_floor_ratio": pooled[
                "exact_to_mapping_floor_ratio"
            ],
            "required_ratio": pooled["required_ratio"],
            "passed": pooled["materiality_gate_passed"],
        },
        "formal_decision": (
            "interaction_not_resolvable_above_numerical_floor_no_refit"
        ),
        "model_refit_executed": False,
        "validation_response_read": False,
        "sealed_test_response_read": False,
        "scientific_interpretation": [
            (
                "The acquisition succeeded algebraically: four complete "
                "fidelity lattices span affine geometry rank 3."
            ),
            (
                "The primary physical interaction did not clear the "
                "predeclared numerical materiality floor on either primary "
                "geometry or after pooling."
            ),
            (
                "Therefore contrast identifiability is necessary but not "
                "sufficient for physically resolvable attribution in this "
                "CFD family."
            ),
            (
                "This is a prospectively obtained negative result, not "
                "evidence that a different neural architecture would solve "
                "the problem."
            ),
        ],
        "claim_scope": {
            "unsupported_claim": (
                "a physically resolved prediction-attribution gap"
            ),
            "supported_interpretation": (
                "algebraic identifiability must be followed by a numerical "
                "resolvability gate before interaction attribution"
            ),
        },
        "next_research_priority": [
            (
                "Do not spend effort changing the neural model while the "
                "target contrast remains below the numerical evidence floor."
            ),
            (
                "First test the identifiability-resolvability hierarchy on "
                "additional already-computed complete lattices and one "
                "independent public CFD family."
            ),
            (
                "Only if at least one interaction clears a prospectively "
                "frozen materiality gate should a new model-attribution "
                "comparison be run."
            ),
        ],
        "input_sha256": {
            path.name: sha256(path)
            for path in (INTEGRATION, FREEZE, IMPLEMENTATION, INPUTS, MATERIALITY)
        },
    }
    OUT_JSON.write_text(json.dumps(decision, indent=2) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Protocol 179 scientific decision\n\n"
        "## Outcome\n\n"
        "The second-stage computation was technically and scientifically "
        "successful: all three new labels passed the unchanged steady-response "
        "platform, producing four strict complete lattices with affine geometry "
        "rank 3. The subsequent pre-fit materiality gate nevertheless failed.\n\n"
        f"- Pooled exact mesh-closure interaction RMS: "
        f"`{pooled['exact_velocity_interaction_rms']:.10g}`.\n"
        f"- Pooled frozen mapping-floor RMS: "
        f"`{pooled['mapping_floor_velocity_rms']:.10g}`.\n"
        f"- Signal/floor ratio: "
        f"`{pooled['exact_to_mapping_floor_ratio']:.6f}`.\n"
        f"- Predeclared requirement: `>= {pooled['required_ratio']:.1f}`.\n\n"
        "The fixed A-C model refit was therefore not run. No validation response "
        "and no sealed test response was opened.\n\n"
        "## Scientific meaning\n\n"
        "The positive prediction-attribution-gap claim is not supported. The "
        "defensible finding is narrower but useful: completing the factorial "
        "contrast support repaired algebraic identifiability, yet the resulting "
        "physical interaction remained below the predeclared numerical "
        "resolvability threshold. In this CFD family, identifiability is "
        "necessary but not sufficient for trustworthy attribution.\n\n"
        "## Recommendation\n\n"
        "Do not tune or replace the neural network to rescue this endpoint. "
        "First establish the proposed two-gate hierarchy—contrast "
        "identifiability, then numerical resolvability—on additional "
        "already-computed complete lattices and an independent public CFD "
        "family. A new attribution fit is justified only after a prospectively "
        "defined target clears the numerical materiality gate.\n",
        encoding="utf-8",
    )
    print(json.dumps(decision, indent=2))


if __name__ == "__main__":
    main()
