#!/usr/bin/env python3
"""Freeze operational details of Protocol-179 interaction evaluation."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PARENT = ROOT / "184_frozen_protocol179_model_comparison.json"
MAPPING = ROOT / "114_geometry_level_mapping_operators.csv"
METRIC_CORE = ROOT / "53_validate_field_metric_pipeline.py"
OUT = ROOT / "184a_frozen_interaction_evaluation_implementation.json"
HASH = ROOT / "184a_interaction_evaluation_sha256.json"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def main() -> None:
    payload = {
        "protocol_id": "184a-v1-20260728",
        "status": "frozen_before_model_refit_or_validation_response_access",
        "parent_protocol_sha256": sha256(PARENT),
        "primary_geometries": [
            "alph15-7929-2024",
            "alph15-13929-2024",
        ],
        "common_support": {
            "source_fidelity": (
                "the source_fidelity_id already frozen per geometry in "
                "114_geometry_level_mapping_operators.csv; F100 for both "
                "primary geometries"
            ),
            "eligible_nodes": (
                "reuse the frozen DNS-to-CFD validity mask from "
                "53_validate_field_metric_pipeline.py"
            ),
            "weights": "physical cell volumes of the source support",
            "mapping": (
                "SciPy LinearNDInterpolator in physical x-y coordinates, with "
                "the same periodic left/right copies used by "
                "mapping_round_trip; intersect finite support across all four "
                "fidelities"
            ),
        },
        "exact_contrast": {
            "velocity": "F110 - F100 - F010 + F000 on common support",
            "pressure": (
                "volume-gauge each mapped pressure field on the common finite "
                "support, then F110 - F100 - F010 + F000"
            ),
            "reattachment": "x_r(F110)-x_r(F100)-x_r(F010)+x_r(F000)",
        },
        "model_interaction": {
            "field": (
                "the interaction head evaluated directly on the frozen common "
                "support and transformed back with training target scale; no "
                "base mean is added"
            ),
            "reattachment": (
                "the interaction logit head is descriptive only because the "
                "sigmoid aggregate is nonlinear"
            ),
        },
        "metrics": {
            "velocity_interaction_nrmse": (
                "sqrt(sum(V*|pred-exact|^2)/sum(V*|exact|^2)), pooled over "
                "the two primary geometries"
            ),
            "amplitude_ratio": (
                "pooled volume-weighted RMS(pred)/RMS(exact)"
            ),
            "spatial_correlation": (
                "pooled volume-weighted centered vector correlation"
            ),
            "pressure_interaction_nrmse": (
                "secondary pooled volume-weighted gauge-free NRMSE"
            ),
        },
        "materiality_floor": {
            "archived_field_floor_per_geometry": (
                "archived mapping_round_trip_velocity_relative_L2 multiplied "
                "by the volume-weighted RMS magnitude of the F100 support "
                "velocity"
            ),
            "pooled_floor": (
                "root-volume-weighted pooling of the two per-geometry floors"
            ),
            "gate": "pooled exact-interaction RMS >= 3*pooled mapping floor",
            "repeatability_note": (
                "no independent full-field repeatability pair was archived; "
                "no unregistered repeatability estimate is introduced"
            ),
        },
        "prediction_evaluation": {
            "validation": "unchanged four development DNS cases",
            "checkpoint": "lowest frozen primary validation score",
            "matched_pairing": "same seed and same K_DNS=12 reveal order",
            "test_access": "prohibited",
        },
        "input_sha256": {
            PARENT.name: sha256(PARENT),
            MAPPING.name: sha256(MAPPING),
            METRIC_CORE.name: sha256(METRIC_CORE),
        },
        "validation_response_read": False,
        "test_response_read": False,
    }
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    HASH.write_text(
        json.dumps({OUT.name: sha256(OUT)}, indent=2) + "\n",
        encoding="utf-8",
    )
    print(HASH.read_text(encoding="utf-8"), end="")


if __name__ == "__main__":
    main()
