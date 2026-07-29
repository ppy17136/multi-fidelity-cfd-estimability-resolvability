#!/usr/bin/env python3
"""Post-Protocol exploratory interaction-materiality audit on four lattices.

This analysis does not alter the frozen Protocol-179 primary endpoint.  It
extends the same exact 2x2 mixed contrast to every strict complete lattice
identified before the refit decision and reports both archived and conservative
mapping floors for the one geometry whose comparison support was refined.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
RUNNER = ROOT / "186_run_protocol179_fixed_model_comparison.py"
INDEX = ROOT / "185_protocol179_cfd_source_cache_index.csv"
DNS_INPUT = ROOT / "109_dns_metric_input_manifest.csv"
OPERATORS = ROOT / "114_geometry_level_mapping_operators.csv"
REFINED = ROOT / "113_refined_mapping_support_metrics.csv"
INTEGRATION = ROOT / "183_second_stage_integration_audit.json"
QOI = ROOT / "185_protocol179_resolved_lower_wall_qoi.csv"
OUT_CSV = ROOT / "188_exploratory_four_lattice_materiality.csv"
OUT_JSON = ROOT / "188_exploratory_four_lattice_materiality.json"
OUT_MD = ROOT / "188_exploratory_four_lattice_materiality.md"
MARKER = ROOT / "EXPLORATORY_FOUR_LATTICE_MATERIALITY_188_COMPLETE.ok"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


R = load_module("runner186_for_audit188", RUNNER)
CORE = R.CORE


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def weighted_rms(value: np.ndarray, weight: np.ndarray) -> float:
    squared = np.sum(value**2, axis=1)
    return float(np.sqrt(np.sum(weight * squared) / np.sum(weight)))


def pooled(rows: list[dict], floor_key: str) -> dict:
    volume = sum(float(row["support_volume"]) for row in rows)
    exact_energy = sum(
        float(row["support_volume"])
        * float(row["exact_velocity_interaction_rms"]) ** 2
        for row in rows
    )
    floor_energy = sum(
        float(row["support_volume"]) * float(row[floor_key]) ** 2
        for row in rows
    )
    exact = math.sqrt(exact_energy / volume)
    floor = math.sqrt(floor_energy / volume)
    ratio = exact / floor
    return {
        "geometry_count": len(rows),
        "support_volume": volume,
        "exact_velocity_interaction_rms": exact,
        "mapping_floor_velocity_rms": floor,
        "exact_to_mapping_floor_ratio": ratio,
        "required_ratio": 3.0,
        "materiality_gate_passed": ratio >= 3.0,
    }


def main() -> None:
    integration = json.loads(INTEGRATION.read_text(encoding="utf-8"))
    cases = list(integration["strict_complete_lattices"])
    if len(cases) != 4:
        raise RuntimeError(f"expected four strict complete lattices, got {cases}")

    cache = {
        (row["case_id"], row["fidelity_id"]): row
        for row in read_csv(INDEX)
    }
    dns = {
        (row["case_id"], row["fidelity_id"]): row
        for row in read_csv(DNS_INPUT)
    }
    operators = {row["case_id"]: row for row in read_csv(OPERATORS)}
    refined = {
        row["case_id"]: row for row in read_csv(REFINED)
    }
    qoi = {
        (row["case_id"], row["fidelity_id"]): float(
            row["CFD_reattachment_x_over_H"]
        )
        for row in read_csv(QOI)
        if row["fidelity_family"] == "steady_RANS"
    }

    rows: list[dict] = []
    for case_id in cases:
        operator = operators[case_id]
        source_fid = operator["source_fidelity_id"]
        source = cache[(case_id, source_fid)]
        source_dir = Path(source["source_metric_directory"])
        centres = CORE.read_openfoam_internal(source_dir / "C", 3)
        volumes = CORE.read_openfoam_internal(source_dir / "Vc", 1)
        source_velocity = CORE.read_openfoam_internal(source_dir / "U", 3)[
            :, :2
        ]

        metric = dns.get((case_id, source_fid)) or dns[(case_id, "F000")]
        x_dns, y_dns, dns_field = CORE.read_dns(
            Path(metric["dns_extracted_path"])
        )
        cache_data = np.load(source["cache_path"])
        alpha = float(cache_data["geometry"][0])
        lx = float(cache_data["Lx_over_H"])
        _, eligible = CORE.interpolate_dns_to_points(
            x_dns,
            y_dns,
            dns_field,
            centres,
            [2, 3, 4, 5],
            alpha,
            lx,
        )
        physical_support_count = int(np.sum(eligible))
        query = centres[eligible]
        weight = volumes[eligible]

        mapped: dict[str, np.ndarray] = {}
        for fidelity_id in ("F000", "F100", "F010", "F110"):
            local = cache[(case_id, fidelity_id)]
            directory = Path(local["source_metric_directory"])
            local_centres = CORE.read_openfoam_internal(directory / "C", 3)
            velocity = CORE.read_openfoam_internal(directory / "U", 3)[:, :2]
            pressure = CORE.read_openfoam_internal(directory / "p", 1)
            mapped[fidelity_id] = R.periodic_linear_map(
                local_centres,
                np.column_stack((velocity, pressure)),
                query,
                lx,
            )

        finite = np.ones(len(query), dtype=bool)
        for value in mapped.values():
            finite &= np.isfinite(value).all(axis=1)
        query = query[finite]
        weight = weight[finite]
        for fidelity_id in mapped:
            mapped[fidelity_id] = mapped[fidelity_id][finite]
            mapped[fidelity_id][:, 2] -= np.average(
                mapped[fidelity_id][:, 2], weights=weight
            )

        exact = (
            mapped["F110"]
            - mapped["F100"]
            - mapped["F010"]
            + mapped["F000"]
        )
        exact_rms = weighted_rms(exact[:, :2], weight)
        source_rms = weighted_rms(source_velocity[eligible][finite], weight)
        archived_relative = float(
            operator["mapping_round_trip_velocity_relative_L2"]
        )
        conservative_relative = archived_relative
        refined_row = refined.get(case_id)
        uses_refined_operator = (
            operator["operator_type"]
            == "deterministic_linear_comparison_support_refinement"
        )
        if uses_refined_operator:
            if refined_row is None:
                raise RuntimeError(
                    f"missing pre-refinement audit for {case_id}"
                )
            conservative_relative = max(
                archived_relative,
                float(refined_row["original_mapping_error"]),
            )
        archived_floor = archived_relative * source_rms
        conservative_floor = conservative_relative * source_rms
        exact_qoi = (
            qoi[(case_id, "F110")]
            - qoi[(case_id, "F100")]
            - qoi[(case_id, "F010")]
            + qoi[(case_id, "F000")]
        )
        rows.append(
            {
                "case_id": case_id,
                "source_fidelity_id": source_fid,
                "operator_type": operator["operator_type"],
                "operator_support_point_count": int(
                    operator["support_point_count"]
                ),
                "physical_cell_support_points_before_finite_intersection": (
                    physical_support_count
                ),
                "finite_common_support_points": int(np.sum(finite)),
                "finite_support_fraction": float(np.mean(finite)),
                "support_volume": float(np.sum(weight)),
                "exact_velocity_interaction_rms": exact_rms,
                "source_velocity_rms": source_rms,
                "archived_mapping_relative_L2": archived_relative,
                "archived_mapping_floor_velocity_rms": archived_floor,
                "exact_to_archived_floor_ratio": exact_rms / archived_floor,
                "archived_three_times_floor_passed": (
                    exact_rms >= 3.0 * archived_floor
                ),
                "conservative_mapping_relative_L2": conservative_relative,
                "conservative_mapping_floor_velocity_rms": conservative_floor,
                "exact_to_conservative_floor_ratio": (
                    exact_rms / conservative_floor
                ),
                "conservative_three_times_floor_passed": (
                    exact_rms >= 3.0 * conservative_floor
                ),
                "exact_pressure_interaction_rms_gauge_free": weighted_rms(
                    exact[:, 2:3], weight
                ),
                "exact_reattachment_mixed_contrast_over_H": exact_qoi,
                "support_note": (
                    "physical source-cell centres and volumes; archived "
                    "refined comparison error is reported separately from "
                    "the conservative pre-refinement error"
                    if uses_refined_operator
                    else "frozen highest-resolution physical-cell support"
                ),
            }
        )

    with OUT_CSV.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    original_cases = {"alph15-7929-2024", "alph15-13929-2024"}
    original = [row for row in rows if row["case_id"] in original_cases]
    completion = [row for row in rows if row["case_id"] not in original_cases]
    payload = {
        "status": "post_protocol_exploratory_not_a_replacement_primary_endpoint",
        "strict_complete_lattices": cases,
        "per_geometry": rows,
        "pooled_all_four_archived_floor": pooled(
            rows, "archived_mapping_floor_velocity_rms"
        ),
        "pooled_all_four_conservative_floor": pooled(
            rows, "conservative_mapping_floor_velocity_rms"
        ),
        "pooled_original_two_archived_floor": pooled(
            original, "archived_mapping_floor_velocity_rms"
        ),
        "pooled_completion_two_archived_floor": pooled(
            completion, "archived_mapping_floor_velocity_rms"
        ),
        "interpretation_rule": (
            "A ratio below 3 means exact factorial identifiability has been "
            "restored but the physical interaction remains unresolved above "
            "the archived numerical mapping evidence. A ratio at or above 3 "
            "is exploratory evidence only and requires independent numerical "
            "repeatability or refinement confirmation before attribution."
        ),
        "validation_response_read": False,
        "test_response_read": False,
        "input_sha256": {
            path.name: sha256(path)
            for path in (
                INDEX,
                DNS_INPUT,
                OPERATORS,
                REFINED,
                INTEGRATION,
                QOI,
            )
        },
    }
    OUT_JSON.write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    all_archived = payload["pooled_all_four_archived_floor"]
    all_conservative = payload["pooled_all_four_conservative_floor"]
    lines = [
        "# Exploratory four-lattice interaction-materiality audit",
        "",
        "This post-Protocol analysis does not replace the frozen Protocol-179 "
        "primary decision and opens no validation or sealed-test response.",
        "",
        "| Geometry | Exact RMS | Archived floor | Ratio | Conservative ratio |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['case_id']} | "
            f"{row['exact_velocity_interaction_rms']:.8g} | "
            f"{row['archived_mapping_floor_velocity_rms']:.8g} | "
            f"{row['exact_to_archived_floor_ratio']:.3f} | "
            f"{row['exact_to_conservative_floor_ratio']:.3f} |"
        )
    lines.extend(
        [
            "",
            f"- Four-lattice pooled archived-floor ratio: "
            f"`{all_archived['exact_to_mapping_floor_ratio']:.6f}`.",
            f"- Four-lattice pooled conservative-floor ratio: "
            f"`{all_conservative['exact_to_mapping_floor_ratio']:.6f}`.",
            "- Required physical-materiality ratio: `>= 3`.",
            "- Exact contrast estimability and numerical resolvability are "
            "reported as separate gates.",
            "",
        ]
    )
    OUT_MD.write_text("\n".join(lines), encoding="utf-8")
    MARKER.write_text(
        "status=complete_post_protocol_exploratory\n"
        f"all_four_archived_ratio="
        f"{all_archived['exact_to_mapping_floor_ratio']:.12g}\n"
        f"all_four_conservative_ratio="
        f"{all_conservative['exact_to_mapping_floor_ratio']:.12g}\n"
        "validation response read: no\n"
        "test response read: no\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

