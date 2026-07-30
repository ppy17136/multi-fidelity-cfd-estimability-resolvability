"""Contrast-level mapping-pathway sensitivity for two frozen CFD geometries.

Pathway A reproduces the archived periodic linear-Delaunay mapping. Pathway B
uses a separately implemented periodic k-nearest-neighbour inverse-distance
map in transformed coordinates. The cell-wise pathway differences are
propagated to the declared mixed contrast with the triangle inequality.
"""

from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
from scipy.interpolate import LinearNDInterpolator
from scipy.spatial import cKDTree


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
INDEX = ROOT / "185_protocol179_cfd_source_cache_index.csv"
DNS_INPUT = ROOT / "109_dns_metric_input_manifest.csv"
OPERATORS = ROOT / "114_geometry_level_mapping_operators.csv"
ARCHIVED_FOUR = ROOT / "188_exploratory_four_lattice_materiality.csv"
CORE_SOURCE = ROOT / "53_validate_field_metric_pipeline.py"
PROTOCOL = HERE / "05_FROZEN_REAL_CFD_MAPPING_REPAIR_IMPLEMENTATION.md"
PRIMARY_GEOMETRIES = ("alph15-7929-2024", "alph15-13929-2024")
FIDELITIES = ("F000", "F100", "F010", "F110")
COEFFICIENTS = {
    "F000": 1.0,
    "F100": -1.0,
    "F010": -1.0,
    "F110": 1.0,
}
PRIMARY_K = 8
ROBUSTNESS_K = (4, 8, 16)
IDW_POWER = 2.0


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


CORE = load_module("metric_core_v9_mapping", CORE_SOURCE)


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
    squared = np.sum(np.asarray(value, dtype=float) ** 2, axis=1)
    return float(np.sqrt(np.sum(weight * squared) / np.sum(weight)))


def periodic_linear_map(
    centres: np.ndarray,
    values: np.ndarray,
    query: np.ndarray,
    lx: float,
) -> np.ndarray:
    """Archived pathway A, reproduced from Protocol-179."""
    unique_x = len(np.unique(np.round(centres[:, 0], 12)))
    strip = 2.1 * lx / max(2, unique_x)
    left = centres[:, 0] < strip
    right = centres[:, 0] > lx - strip
    points = np.vstack(
        (
            centres[:, :2],
            centres[left, :2] + np.asarray([lx, 0.0]),
            centres[right, :2] - np.asarray([lx, 0.0]),
        )
    )
    mapped_values = np.concatenate((values, values[left], values[right]), axis=0)
    interpolator = LinearNDInterpolator(points, mapped_values, fill_value=np.nan)
    return np.asarray(interpolator(query[:, :2]))


def periodic_idw_map(
    coordinates: np.ndarray,
    values: np.ndarray,
    query: np.ndarray,
    *,
    neighbours: int,
    power: float,
) -> np.ndarray:
    """Independent pathway B in periodic transformed coordinates."""
    coordinates = np.asarray(coordinates, dtype=float)
    values = np.asarray(values, dtype=float)
    query = np.asarray(query, dtype=float)
    shifts = (-1.0, 0.0, 1.0)
    points = np.vstack(
        [coordinates + np.asarray([shift, 0.0]) for shift in shifts]
    )
    repeated_values = np.vstack([values for _ in shifts])
    tree = cKDTree(points)
    distances, indices = tree.query(
        query,
        k=min(neighbours, len(points)),
        workers=-1,
    )
    if distances.ndim == 1:
        distances = distances[:, None]
        indices = indices[:, None]
    selected = repeated_values[indices]
    output = np.empty((len(query), values.shape[1]), dtype=float)
    exact = distances <= 1.0e-12
    exact_rows = np.any(exact, axis=1)
    for row in np.flatnonzero(exact_rows):
        output[row] = selected[row, exact[row]].mean(axis=0)
    nonexact_rows = np.flatnonzero(~exact_rows)
    if len(nonexact_rows):
        local_distance = distances[nonexact_rows]
        weights = local_distance ** (-power)
        weights /= weights.sum(axis=1, keepdims=True)
        output[nonexact_rows] = np.einsum(
            "nk,nkd->nd", weights, selected[nonexact_rows]
        )
    return output


def evaluate_geometry(
    case_id: str,
    neighbours: int,
    cache: dict[tuple[str, str], dict[str, str]],
    dns: dict[tuple[str, str], dict[str, str]],
    operators: dict[str, dict[str, str]],
    archived: dict[str, dict[str, str]],
) -> tuple[dict, list[dict]]:
    operator = operators[case_id]
    if operator["source_fidelity_id"] != "F100":
        raise RuntimeError(f"unexpected common support for {case_id}")
    source = cache[(case_id, "F100")]
    source_dir = Path(source["source_metric_directory"])
    source_cache = np.load(source["cache_path"])
    centres = CORE.read_openfoam_internal(source_dir / "C", 3)
    volumes = CORE.read_openfoam_internal(source_dir / "Vc", 1)
    if len(centres) != len(source_cache["coord"]):
        raise RuntimeError(f"source coordinate count mismatch for {case_id}")

    metric = dns[(case_id, "F100")]
    x_dns, y_dns, dns_field = CORE.read_dns(Path(metric["dns_extracted_path"]))
    _, eligible = CORE.interpolate_dns_to_points(
        x_dns,
        y_dns,
        dns_field,
        centres,
        [2, 3, 4, 5],
        float(source_cache["geometry"][0]),
        float(source_cache["Lx_over_H"]),
    )
    if int(np.sum(eligible)) != int(operator["support_point_count"]):
        raise RuntimeError(f"frozen support mismatch for {case_id}")

    query_physical = centres[eligible]
    query_transformed = np.asarray(source_cache["coord"][eligible], dtype=float)
    weight = np.asarray(volumes[eligible], dtype=float)
    lx = float(source_cache["Lx_over_H"])

    pathway_a: dict[str, np.ndarray] = {}
    pathway_b: dict[str, np.ndarray] = {}
    for fidelity_id in FIDELITIES:
        row = cache[(case_id, fidelity_id)]
        directory = Path(row["source_metric_directory"])
        local_centres = CORE.read_openfoam_internal(directory / "C", 3)
        local_velocity = CORE.read_openfoam_internal(directory / "U", 3)[:, :2]
        local_cache = np.load(row["cache_path"])
        cached_velocity = np.asarray(local_cache["target"][:, :2], dtype=float)
        if not np.allclose(local_velocity, cached_velocity, rtol=2e-6, atol=2e-7):
            raise RuntimeError(f"cache/source velocity mismatch: {case_id}/{fidelity_id}")
        pathway_a[fidelity_id] = periodic_linear_map(
            local_centres, local_velocity, query_physical, lx
        )
        pathway_b[fidelity_id] = periodic_idw_map(
            np.asarray(local_cache["coord"], dtype=float),
            cached_velocity,
            query_transformed,
            neighbours=neighbours,
            power=IDW_POWER,
        )

    finite = np.ones(len(query_physical), dtype=bool)
    for field in (*pathway_a.values(), *pathway_b.values()):
        finite &= np.isfinite(field).all(axis=1)
    finite &= np.isfinite(weight) & (weight > 0)
    if not np.any(finite):
        raise RuntimeError(f"no finite common support for {case_id}")
    weight = weight[finite]
    for fidelity_id in FIDELITIES:
        pathway_a[fidelity_id] = pathway_a[fidelity_id][finite]
        pathway_b[fidelity_id] = pathway_b[fidelity_id][finite]

    interaction_a = sum(
        COEFFICIENTS[fid] * pathway_a[fid] for fid in FIDELITIES
    )
    interaction_b = sum(
        COEFFICIENTS[fid] * pathway_b[fid] for fid in FIDELITIES
    )
    cell_rows = []
    cell_bounds = []
    for fidelity_id in FIDELITIES:
        sensitivity = weighted_rms(
            pathway_a[fidelity_id] - pathway_b[fidelity_id], weight
        )
        cell_bounds.append(abs(COEFFICIENTS[fidelity_id]) * sensitivity)
        cell_rows.append(
            {
                "case_id": case_id,
                "neighbours": neighbours,
                "fidelity_id": fidelity_id,
                "contrast_coefficient": COEFFICIENTS[fidelity_id],
                "finite_common_support_points": int(np.sum(finite)),
                "cell_pathway_sensitivity_velocity_rms": sensitivity,
                "cell_pathway_a_velocity_rms": weighted_rms(
                    pathway_a[fidelity_id], weight
                ),
                "cell_pathway_b_velocity_rms": weighted_rms(
                    pathway_b[fidelity_id], weight
                ),
            }
        )

    propagated_bound = float(sum(cell_bounds))
    observed_pathway_difference = weighted_rms(interaction_a - interaction_b, weight)
    interaction_a_rms = weighted_rms(interaction_a, weight)
    interaction_b_rms = weighted_rms(interaction_b, weight)
    archived_row = archived[case_id]
    result = {
        "case_id": case_id,
        "neighbours": neighbours,
        "primary_setting": neighbours == PRIMARY_K,
        "operator_support_point_count": int(operator["support_point_count"]),
        "finite_common_support_points": int(np.sum(finite)),
        "finite_support_fraction": float(np.mean(finite)),
        "support_volume": float(np.sum(weight)),
        "pathway_a_interaction_velocity_rms": interaction_a_rms,
        "pathway_b_interaction_velocity_rms": interaction_b_rms,
        "observed_interaction_pathway_difference_rms": observed_pathway_difference,
        "descriptive_interaction_to_observed_pathway_difference_ratio": (
            interaction_a_rms / observed_pathway_difference
        ),
        "propagated_cellwise_pathway_sensitivity_bound": propagated_bound,
        "triangle_inequality_verified": observed_pathway_difference <= (
            propagated_bound + 1.0e-14
        ),
        "pathway_a_interaction_to_bound_ratio": (
            interaction_a_rms / propagated_bound
        ),
        "three_times_pathway_bound_passed": (
            interaction_a_rms >= 3.0 * propagated_bound
        ),
        "archived_single_floor_interaction_rms": float(
            archived_row["exact_velocity_interaction_rms"]
        ),
        "archived_single_floor_ratio": float(
            archived_row["exact_to_archived_floor_ratio"]
        ),
        "claim_boundary": (
            "The propagated quantity bounds disagreement between two mapping "
            "implementations; it is not a total numerical-error bound. The "
            "direct pathway-difference ratio is descriptive repeatability "
            "evidence and is not used as a deterministic resolution bound."
        ),
    }
    return result, cell_rows


def main() -> None:
    cache = {
        (row["case_id"], row["fidelity_id"]): row for row in read_csv(INDEX)
    }
    dns = {
        (row["case_id"], row["fidelity_id"]): row for row in read_csv(DNS_INPUT)
    }
    operators = {row["case_id"]: row for row in read_csv(OPERATORS)}
    archived = {row["case_id"]: row for row in read_csv(ARCHIVED_FOUR)}

    geometry_rows = []
    cell_rows = []
    for neighbours in ROBUSTNESS_K:
        for case_id in PRIMARY_GEOMETRIES:
            geometry, cells = evaluate_geometry(
                case_id, neighbours, cache, dns, operators, archived
            )
            geometry_rows.append(geometry)
            cell_rows.extend(cells)

    with (HERE / "05_real_cfd_contrast_mapping_repair.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(geometry_rows[0]))
        writer.writeheader()
        writer.writerows(geometry_rows)
    with (HERE / "05_real_cfd_cellwise_mapping_sensitivity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(cell_rows[0]))
        writer.writeheader()
        writer.writerows(cell_rows)

    primary = [row for row in geometry_rows if row["primary_setting"]]
    payload = {
        "protocol": "V9-20260730",
        "analysis_status": "retrospective archived-field repair",
        "primary_idw_neighbours": PRIMARY_K,
        "robustness_idw_neighbours": list(ROBUSTNESS_K),
        "idw_power": IDW_POWER,
        "primary_geometries": list(PRIMARY_GEOMETRIES),
        "primary_results": primary,
        "all_results": geometry_rows,
        "all_triangle_inequality_checks_passed": all(
            row["triangle_inequality_verified"] for row in geometry_rows
        ),
        "primary_three_times_bound_pass_count": sum(
            row["three_times_pathway_bound_passed"] for row in primary
        ),
        "input_sha256": {
            path.name: sha256(path)
            for path in (
                INDEX,
                DNS_INPUT,
                OPERATORS,
                ARCHIVED_FOUR,
                CORE_SOURCE,
                PROTOCOL,
            )
        },
        "claim_boundary": (
            "This repair makes the cell-to-contrast propagation explicit and "
            "tests an independent mapping pathway. It does not estimate total "
            "CFD, turbulence-model, or physical model discrepancy."
        ),
    }
    (HERE / "05_real_cfd_contrast_mapping_repair.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(payload, indent=2))


if __name__ == "__main__":
    main()

