"""Compare axis-aligned and joint mapping-ensemble uncertainty sets."""

from __future__ import annotations

import csv
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
BASE_IMPL = HERE / "05_run_real_cfd_contrast_mapping_repair.py"
INDEX = ROOT / "185_protocol179_cfd_source_cache_index.csv"
DNS_INPUT = ROOT / "109_dns_metric_input_manifest.csv"
OPERATORS = ROOT / "114_geometry_level_mapping_operators.csv"
PRIMARY_GEOMETRIES = ("alph15-7929-2024", "alph15-13929-2024")
FIDELITIES = ("F000", "F100", "F010", "F110")
COEFFICIENTS = {
    "F000": 1.0,
    "F100": -1.0,
    "F010": -1.0,
    "F110": 1.0,
}
PATHWAYS = ("linear", "idw4", "idw8", "idw16")


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


B = load_module("mapping_repair_base07", BASE_IMPL)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        return list(csv.DictReader(handle))


def prepare_support(
    case_id: str,
    source_fidelity: str,
    cache: dict[tuple[str, str], dict[str, str]],
    dns: dict[tuple[str, str], dict[str, str]],
) -> tuple[np.ndarray, np.ndarray, np.ndarray, float]:
    source = cache[(case_id, source_fidelity)]
    source_dir = Path(source["source_metric_directory"])
    source_cache = np.load(source["cache_path"])
    centres = B.CORE.read_openfoam_internal(source_dir / "C", 3)
    volumes = B.CORE.read_openfoam_internal(source_dir / "Vc", 1)
    metric = dns[(case_id, source_fidelity)]
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
    return (
        centres[eligible],
        np.asarray(source_cache["coord"][eligible], dtype=float),
        np.asarray(volumes[eligible], dtype=float),
        float(source_cache["Lx_over_H"]),
    )


def evaluate_support(
    case_id: str,
    support_name: str,
    source_fidelity: str,
    cache: dict[tuple[str, str], dict[str, str]],
    dns: dict[tuple[str, str], dict[str, str]],
) -> tuple[dict, list[dict]]:
    query_physical, query_transformed, weights, lx = prepare_support(
        case_id, source_fidelity, cache, dns
    )
    fields: dict[str, dict[str, np.ndarray]] = {
        pathway: {} for pathway in PATHWAYS
    }
    for fidelity_id in FIDELITIES:
        row = cache[(case_id, fidelity_id)]
        directory = Path(row["source_metric_directory"])
        centres = B.CORE.read_openfoam_internal(directory / "C", 3)
        velocity = B.CORE.read_openfoam_internal(directory / "U", 3)[:, :2]
        local_cache = np.load(row["cache_path"])
        fields["linear"][fidelity_id] = B.periodic_linear_map(
            centres, velocity, query_physical, lx
        )
        for neighbours in (4, 8, 16):
            fields[f"idw{neighbours}"][fidelity_id] = B.periodic_idw_map(
                np.asarray(local_cache["coord"], dtype=float),
                np.asarray(local_cache["target"][:, :2], dtype=float),
                query_transformed,
                neighbours=neighbours,
                power=2.0,
            )

    finite = np.ones(len(query_physical), dtype=bool)
    for pathway in PATHWAYS:
        for fidelity_id in FIDELITIES:
            finite &= np.isfinite(fields[pathway][fidelity_id]).all(axis=1)
    finite &= np.isfinite(weights) & (weights > 0)
    weights = weights[finite]
    for pathway in PATHWAYS:
        for fidelity_id in FIDELITIES:
            fields[pathway][fidelity_id] = fields[pathway][fidelity_id][finite]

    interactions = {
        pathway: sum(
            COEFFICIENTS[fidelity_id] * fields[pathway][fidelity_id]
            for fidelity_id in FIDELITIES
        )
        for pathway in PATHWAYS
    }
    signal_rms = {
        pathway: B.weighted_rms(interactions[pathway], weights)
        for pathway in PATHWAYS
    }
    interaction_pairs = []
    for first, second in itertools.combinations(PATHWAYS, 2):
        difference = B.weighted_rms(
            interactions[first] - interactions[second], weights
        )
        interaction_pairs.append(
            {
                "case_id": case_id,
                "support": support_name,
                "quantity": "mixed_contrast",
                "first_pathway": first,
                "second_pathway": second,
                "difference_velocity_rms": difference,
            }
        )
    joint_diameter_row = max(
        interaction_pairs, key=lambda row: row["difference_velocity_rms"]
    )
    joint_diameter = float(joint_diameter_row["difference_velocity_rms"])

    cell_diameters = []
    pair_rows = list(interaction_pairs)
    for fidelity_id in FIDELITIES:
        local_pairs = []
        for first, second in itertools.combinations(PATHWAYS, 2):
            difference = B.weighted_rms(
                fields[first][fidelity_id] - fields[second][fidelity_id],
                weights,
            )
            row = {
                "case_id": case_id,
                "support": support_name,
                "quantity": fidelity_id,
                "first_pathway": first,
                "second_pathway": second,
                "difference_velocity_rms": difference,
            }
            pair_rows.append(row)
            local_pairs.append(row)
        cell_diameters.append(
            max(row["difference_velocity_rms"] for row in local_pairs)
        )
    axis_aligned_bound = float(sum(cell_diameters))
    conservative_signal = float(min(signal_rms.values()))
    result = {
        "case_id": case_id,
        "support": support_name,
        "source_fidelity": source_fidelity,
        "finite_common_support_points": int(np.sum(finite)),
        "support_volume": float(np.sum(weights)),
        "linear_interaction_velocity_rms": signal_rms["linear"],
        "minimum_ensemble_interaction_velocity_rms": conservative_signal,
        "maximum_ensemble_interaction_velocity_rms": float(
            max(signal_rms.values())
        ),
        "joint_mapping_ensemble_diameter": joint_diameter,
        "diameter_attaining_pair": (
            f"{joint_diameter_row['first_pathway']}:"
            f"{joint_diameter_row['second_pathway']}"
        ),
        "ensemble_resolution_ratio": conservative_signal / joint_diameter,
        "ensemble_three_times_diameter_passed": (
            conservative_signal >= 3.0 * joint_diameter
        ),
        "axis_aligned_cellwise_ensemble_bound": axis_aligned_bound,
        "axis_aligned_resolution_ratio": (
            conservative_signal / axis_aligned_bound
        ),
        "axis_aligned_three_times_bound_passed": (
            conservative_signal >= 3.0 * axis_aligned_bound
        ),
        "joint_to_axis_aligned_diameter_ratio": (
            joint_diameter / axis_aligned_bound
        ),
        "uncertainty_set_scope": (
            "finite prespecified mapping ensemble versus independent "
            "cell-wise ensemble envelope"
        ),
    }
    return result, pair_rows


def main() -> None:
    cache = {
        (row["case_id"], row["fidelity_id"]): row for row in read_csv(INDEX)
    }
    dns = {
        (row["case_id"], row["fidelity_id"]): row for row in read_csv(DNS_INPUT)
    }
    operators = {row["case_id"]: row for row in read_csv(OPERATORS)}
    rows = []
    pair_rows = []
    for case_id in PRIMARY_GEOMETRIES:
        if operators[case_id]["source_fidelity_id"] != "F100":
            raise RuntimeError(f"unexpected fine support for {case_id}")
        for support_name, source_fidelity in (
            ("fine_F100", "F100"),
            ("shared_F000", "F000"),
        ):
            result, pairs = evaluate_support(
                case_id, support_name, source_fidelity, cache, dns
            )
            rows.append(result)
            pair_rows.extend(pairs)

    with (HERE / "07_joint_mapping_ensemble_audit.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    with (HERE / "07_joint_mapping_ensemble_pairwise.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(pair_rows[0]))
        writer.writeheader()
        writer.writerows(pair_rows)

    summary = {
        "protocol": "V9-20260730",
        "declared_pathways": list(PATHWAYS),
        "results": rows,
        "joint_ensemble_pass_count": sum(
            row["ensemble_three_times_diameter_passed"] for row in rows
        ),
        "axis_aligned_pass_count": sum(
            row["axis_aligned_three_times_bound_passed"] for row in rows
        ),
        "interpretation": (
            "A joint-ensemble pass is robustness only over the prespecified "
            "mapping algorithms. The axis-aligned result remains the more "
            "adversarial independent-cell statement."
        ),
    }
    (HERE / "07_joint_mapping_ensemble_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

