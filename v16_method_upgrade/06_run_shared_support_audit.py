"""Evaluate the four-cell contrast on the shared coarse representable space."""

from __future__ import annotations

import csv
import importlib.util
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE_IMPL = HERE / "05_run_real_cfd_contrast_mapping_repair.py"
INDEX = ROOT / "185_protocol179_cfd_source_cache_index.csv"
DNS_INPUT = ROOT / "109_dns_metric_input_manifest.csv"
FINE_RESULTS = HERE / "05_real_cfd_contrast_mapping_repair.csv"
PROTOCOL = HERE / "06_FROZEN_SHARED_SUPPORT_AUDIT.md"
PRIMARY_GEOMETRIES = ("alph15-7929-2024", "alph15-13929-2024")
FIDELITIES = ("F000", "F100", "F010", "F110")
COEFFICIENTS = {
    "F000": 1.0,
    "F100": -1.0,
    "F010": -1.0,
    "F110": 1.0,
}
K_VALUES = (4, 8, 16)
PRIMARY_K = 8


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_module("mapping_repair_base06", BASE_IMPL)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def evaluate(
    case_id: str,
    neighbours: int,
    cache: dict[tuple[str, str], dict[str, str]],
    dns: dict[tuple[str, str], dict[str, str]],
    fine: dict[str, dict[str, str]],
) -> tuple[dict, list[dict]]:
    source = cache[(case_id, "F000")]
    source_dir = Path(source["source_metric_directory"])
    source_cache = np.load(source["cache_path"])
    centres = B.CORE.read_openfoam_internal(source_dir / "C", 3)
    volumes = B.CORE.read_openfoam_internal(source_dir / "Vc", 1)
    metric = dns[(case_id, "F000")]
    x_dns, y_dns, dns_field = B.CORE.read_dns(Path(metric["dns_extracted_path"]))
    _, eligible = B.CORE.interpolate_dns_to_points(
        x_dns,
        y_dns,
        dns_field,
        centres,
        [2, 3, 4, 5],
        float(source_cache["geometry"][0]),
        float(source_cache["Lx_over_H"]),
    )
    query_physical = centres[eligible]
    query_transformed = np.asarray(source_cache["coord"][eligible], dtype=float)
    weights = np.asarray(volumes[eligible], dtype=float)
    lx = float(source_cache["Lx_over_H"])

    pathway_a = {}
    pathway_b = {}
    for fidelity_id in FIDELITIES:
        row = cache[(case_id, fidelity_id)]
        directory = Path(row["source_metric_directory"])
        local_centres = B.CORE.read_openfoam_internal(directory / "C", 3)
        velocity = B.CORE.read_openfoam_internal(directory / "U", 3)[:, :2]
        local_cache = np.load(row["cache_path"])
        pathway_a[fidelity_id] = B.periodic_linear_map(
            local_centres, velocity, query_physical, lx
        )
        pathway_b[fidelity_id] = B.periodic_idw_map(
            np.asarray(local_cache["coord"], dtype=float),
            np.asarray(local_cache["target"][:, :2], dtype=float),
            query_transformed,
            neighbours=neighbours,
            power=2.0,
        )

    finite = np.ones(len(query_physical), dtype=bool)
    for field in (*pathway_a.values(), *pathway_b.values()):
        finite &= np.isfinite(field).all(axis=1)
    finite &= np.isfinite(weights) & (weights > 0)
    weights = weights[finite]
    for fidelity_id in FIDELITIES:
        pathway_a[fidelity_id] = pathway_a[fidelity_id][finite]
        pathway_b[fidelity_id] = pathway_b[fidelity_id][finite]

    interaction_a = sum(
        COEFFICIENTS[fidelity_id] * pathway_a[fidelity_id]
        for fidelity_id in FIDELITIES
    )
    interaction_b = sum(
        COEFFICIENTS[fidelity_id] * pathway_b[fidelity_id]
        for fidelity_id in FIDELITIES
    )
    cell_rows = []
    cell_sensitivities = []
    for fidelity_id in FIDELITIES:
        sensitivity = B.weighted_rms(
            pathway_a[fidelity_id] - pathway_b[fidelity_id], weights
        )
        cell_sensitivities.append(sensitivity)
        cell_rows.append(
            {
                "case_id": case_id,
                "neighbours": neighbours,
                "fidelity_id": fidelity_id,
                "cell_pathway_sensitivity_velocity_rms": sensitivity,
            }
        )

    interaction_a_rms = B.weighted_rms(interaction_a, weights)
    interaction_b_rms = B.weighted_rms(interaction_b, weights)
    direct_difference = B.weighted_rms(interaction_a - interaction_b, weights)
    propagated_bound = float(sum(cell_sensitivities))
    fine_row = fine[case_id]
    return (
        {
            "case_id": case_id,
            "neighbours": neighbours,
            "primary_setting": neighbours == PRIMARY_K,
            "eligible_coarse_support_points": int(np.sum(eligible)),
            "finite_shared_support_points": int(np.sum(finite)),
            "finite_support_fraction": float(np.mean(finite)),
            "support_volume": float(np.sum(weights)),
            "shared_pathway_a_interaction_velocity_rms": interaction_a_rms,
            "shared_pathway_b_interaction_velocity_rms": interaction_b_rms,
            "shared_observed_pathway_difference_rms": direct_difference,
            "shared_descriptive_interaction_to_pathway_difference_ratio": (
                interaction_a_rms / direct_difference
            ),
            "shared_propagated_cellwise_pathway_bound": propagated_bound,
            "shared_interaction_to_propagated_bound_ratio": (
                interaction_a_rms / propagated_bound
            ),
            "shared_three_times_bound_passed": (
                interaction_a_rms >= 3.0 * propagated_bound
            ),
            "triangle_inequality_verified": direct_difference <= (
                propagated_bound + 1.0e-14
            ),
            "fine_support_interaction_velocity_rms": float(
                fine_row["pathway_a_interaction_velocity_rms"]
            ),
            "shared_to_fine_interaction_rms_ratio": (
                interaction_a_rms
                / float(fine_row["pathway_a_interaction_velocity_rms"])
            ),
            "estimand_scope": "shared coarse representable subspace",
        },
        cell_rows,
    )


def main() -> None:
    cache = {
        (row["case_id"], row["fidelity_id"]): row for row in read_csv(INDEX)
    }
    dns = {
        (row["case_id"], row["fidelity_id"]): row for row in read_csv(DNS_INPUT)
    }
    fine = {
        row["case_id"]: row
        for row in read_csv(FINE_RESULTS)
        if row["primary_setting"].lower() == "true"
    }
    geometry_rows = []
    cell_rows = []
    for neighbours in K_VALUES:
        for case_id in PRIMARY_GEOMETRIES:
            geometry, cells = evaluate(case_id, neighbours, cache, dns, fine)
            geometry_rows.append(geometry)
            cell_rows.extend(cells)

    with (HERE / "06_shared_support_audit.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(geometry_rows[0]))
        writer.writeheader()
        writer.writerows(geometry_rows)
    with (HERE / "06_shared_support_cellwise_sensitivity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cell_rows[0]))
        writer.writeheader()
        writer.writerows(cell_rows)

    primary = [row for row in geometry_rows if row["primary_setting"]]
    summary = {
        "protocol": "V12-20260730", "frozen_parent_protocol": "V9-20260730",
        "primary_neighbours": PRIMARY_K,
        "primary_results": primary,
        "all_results": geometry_rows,
        "primary_pass_count": sum(
            row["shared_three_times_bound_passed"] for row in primary
        ),
        "all_triangle_inequality_checks_passed": all(
            row["triangle_inequality_verified"] for row in geometry_rows
        ),
        "interpretation": (
            "Shared-support and fine-support contrasts are distinct declared "
            "estimands. A shared-support pass cannot validate fine-scale "
            "resolution-exclusive content."
        ),
    }
    (HERE / "06_shared_support_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
