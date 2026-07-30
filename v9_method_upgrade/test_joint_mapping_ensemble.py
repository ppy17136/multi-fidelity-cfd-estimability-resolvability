#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def test_joint_and_axis_results_are_deliberately_distinct():
    data = json.loads(
        (ROOT / "07_joint_mapping_ensemble_audit.json").read_text(
            encoding="utf-8"
        )
    )
    assert data["joint_ensemble_pass_count"] == 4
    assert data["axis_aligned_pass_count"] == 0
    assert len(data["results"]) == 4
    for row in data["results"]:
        assert row["ensemble_resolution_ratio"] >= 3.0
        assert row["axis_aligned_resolution_ratio"] < 3.0
        assert 0.0 < row["joint_to_axis_aligned_diameter_ratio"] < 1.0
        assert row["diameter_attaining_pair"]


def test_finite_ensemble_scope_is_explicit():
    text = (ROOT / "07_FROZEN_JOINT_MAPPING_ENSEMBLE_AUDIT.md").read_text(
        encoding="utf-8"
    )
    lowered = text.lower()
    assert "finite" in lowered
    assert "total cfd" in lowered
