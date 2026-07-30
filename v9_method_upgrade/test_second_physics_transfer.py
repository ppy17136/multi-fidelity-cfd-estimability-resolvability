#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def test_transfer_designs_are_complete_but_gate2_is_not_overclaimed():
    data = json.loads(
        (ROOT / "09_second_physics_transfer_audit.json").read_text(
            encoding="utf-8"
        )
    )
    assert data["conditions"] == 3
    assert data["complete_four_cell_designs"] == 3
    assert data["gate2_assessed"] is False
    assert len(data["results"]) == 3
    assert all(row["gate1_estimable"] for row in data["results"])
    assert all("not assessed" in row["gate2_status"] for row in data["results"])

