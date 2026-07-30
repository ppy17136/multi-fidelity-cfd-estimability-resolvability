#!/usr/bin/env python3
"""Retrospective mesh-by-closure transfer audit for archived pump CFD."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
PROJECT_ROOT = ROOT.parents[1]
CFD = (
    PROJECT_ROOT
    / "fifth_paper_sim_exp_20260619"
    / "cfd_3d_circular"
)

KE_FILE = (
    CFD
    / "method_upgrade_case_disjoint_response_consistent_20260710"
    / "fixedRe100k_patchK_grid_sensitivity_20260710.csv"
)
SST_FILE = (
    CFD
    / "turbulence_sensitivity_sst_fine_Re100k_v141"
    / "fine_sst_audit_v141.csv"
)
PROTOCOL_FILE = ROOT / "09_FROZEN_SECOND_PHYSICS_TRANSFER_AUDIT.md"


CONDITION_MAP = {
    "normal_dp050": "normal",
    "severe_inlet_dp050_e025": "severe_inlet",
    "severe_outlet_dp050_e025": "severe_outlet",
}


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()


def load_ke() -> dict[tuple[str, str], dict]:
    rows = {}
    with KE_FILE.open(newline="", encoding="utf-8-sig") as stream:
        for raw in csv.DictReader(stream):
            if raw["mesh_level"] not in {"medium", "fine"}:
                continue
            condition = CONDITION_MAP[raw["condition"]]
            rows[(condition, raw["mesh_level"])] = {
                "K": float(raw["K_patch"]),
                "cells": int(raw["cells"]),
                "case": raw["case"],
                "ended": True,
                "mass_imbalance_pct": None,
            }
    return rows


def load_sst() -> dict[tuple[str, str], dict]:
    rows = {}
    with SST_FILE.open(newline="", encoding="utf-8-sig") as stream:
        for raw in csv.DictReader(stream):
            if raw["model"] != "kOmegaSST":
                continue
            if raw["mesh_level"] not in {"medium", "fine"}:
                continue
            rows[(raw["condition"], raw["mesh_level"])] = {
                "K": float(raw["K_boundary"]),
                "cells": int(raw["cells"]),
                "case": raw["case"],
                "ended": raw["ended"].strip().lower() == "yes",
                "mass_imbalance_pct": abs(float(raw["mass_imbalance_pct"])),
            }
    return rows


ke = load_ke()
sst = load_sst()
results = []

for condition in ("normal", "severe_inlet", "severe_outlet"):
    cells = {
        ("medium", "kEpsilon"): ke[(condition, "medium")],
        ("fine", "kEpsilon"): ke[(condition, "fine")],
        ("medium", "kOmegaSST"): sst[(condition, "medium")],
        ("fine", "kOmegaSST"): sst[(condition, "fine")],
    }
    assert len(cells) == 4
    assert all(item["ended"] for item in cells.values())
    assert all(
        item["mass_imbalance_pct"] is None
        or item["mass_imbalance_pct"] < 1.0e-4
        for item in cells.values()
    )

    km_ke = cells[("medium", "kEpsilon")]["K"]
    kf_ke = cells[("fine", "kEpsilon")]["K"]
    km_sst = cells[("medium", "kOmegaSST")]["K"]
    kf_sst = cells[("fine", "kOmegaSST")]["K"]
    contrast = (kf_sst - kf_ke) - (km_sst - km_ke)
    mesh_shift_ke = kf_ke - km_ke
    mesh_shift_sst = kf_sst - km_sst
    closure_shift_medium = km_sst - km_ke
    closure_shift_fine = kf_sst - kf_ke
    mean_abs_k = sum(abs(item["K"]) for item in cells.values()) / 4.0

    results.append(
        {
            "condition": condition,
            "K_medium_kEpsilon": km_ke,
            "K_fine_kEpsilon": kf_ke,
            "K_medium_kOmegaSST": km_sst,
            "K_fine_kOmegaSST": kf_sst,
            "mesh_shift_kEpsilon": mesh_shift_ke,
            "mesh_shift_kOmegaSST": mesh_shift_sst,
            "closure_shift_medium": closure_shift_medium,
            "closure_shift_fine": closure_shift_fine,
            "mesh_by_closure_contrast": contrast,
            "absolute_contrast": abs(contrast),
            "contrast_relative_to_four_cell_mean_abs_pct": (
                100.0 * abs(contrast) / mean_abs_k
            ),
            "gate1_estimable": True,
            "gate2_status": "not assessed: no independent contrast-level K error bound",
            "maximum_reported_sst_mass_imbalance_pct": max(
                item["mass_imbalance_pct"] or 0.0 for item in cells.values()
            ),
        }
    )

csv_path = ROOT / "09_second_physics_transfer_audit.csv"
with csv_path.open("w", newline="", encoding="utf-8") as stream:
    writer = csv.DictWriter(stream, fieldnames=list(results[0]))
    writer.writeheader()
    writer.writerows(results)

summary = {
    "protocol": "V9-20260730",
    "application": "obstructed circular internal flow",
    "fidelity_axes": ["mesh: medium/fine", "closure: kEpsilon/kOmegaSST"],
    "response": "boundary pressure-loss coefficient K",
    "conditions": len(results),
    "complete_four_cell_designs": sum(r["gate1_estimable"] for r in results),
    "results": results,
    "gate2_assessed": False,
    "reason_gate2_not_assessed": (
        "No independent contrast-level pressure-loss error bound is available "
        "for all four archived cells."
    ),
    "recommended_placement": (
        "Brief transfer paragraph/table in the main paper; complete audit in "
        "supplementary material."
    ),
    "input_sha256": {
        KE_FILE.name: sha256(KE_FILE),
        SST_FILE.name: sha256(SST_FILE),
        PROTOCOL_FILE.name: sha256(PROTOCOL_FILE),
    },
}

(ROOT / "09_second_physics_transfer_audit.json").write_text(
    json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
    encoding="utf-8",
)
print(json.dumps(summary, indent=2, ensure_ascii=False))

