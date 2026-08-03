"""Symmetric coherent-versus-independent mapping-set comparison.

For each geometry/support audit, this script constructs:

1. a coherent set of four interaction fields, one per prescribed mapping
   pipeline applied consistently to all four fidelity cells; and
2. an independently recombined set of 4^4 = 256 interaction fields, allowing
   each fidelity cell to select any of the same four mapping pipelines.

Both sets are evaluated with exactly the same functionals:

    S(U) = min_{C in U} ||C||_W
    D(U) = max_{C,C' in U} ||C-C'||_W
    R(U) = S(U) / D(U)

The result isolates the effect of cross-cell coherence from the earlier
triangle-envelope looseness and numerator mismatch.
"""

from __future__ import annotations

import csv
import importlib.util
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
BASE = HERE / "07_run_joint_mapping_ensemble_audit.py"
CACHE_DIR = HERE / "derived_mapping_vectors_v12"


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


J = load_module("mapping_symmetric_base07b", BASE)
B = J.B


def weighted_vector(field: np.ndarray, weights: np.ndarray) -> np.ndarray:
    scale = np.sqrt(weights / np.sum(weights))
    return (field * scale[:, None]).reshape(-1)


def set_metrics(vectors: np.ndarray) -> dict:
    """Return S, D, and R from weighted flattened field vectors."""
    x = np.asarray(vectors, dtype=np.float64)
    gram = x @ x.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = (
        norm2[:, None]
        + norm2[None, :]
        - 2.0 * gram
    )
    np.maximum(distance2, 0.0, out=distance2)
    min_index = int(np.argmin(norm2))
    maximum_flat = int(np.argmax(distance2))
    first, second = np.unravel_index(maximum_flat, distance2.shape)
    gram_signal = float(np.sqrt(norm2[min_index]))
    gram_diameter = float(np.sqrt(distance2[first, second]))
    signal = float(np.linalg.norm(x[min_index]))
    diameter = float(np.linalg.norm(x[first] - x[second]))
    return {
        "set_size": int(len(x)),
        "minimum_norm": signal,
        "diameter": diameter,
        "ratio": signal / diameter if diameter > 0 else float("inf"),
        "three_times_diameter_passed": bool(signal >= 3.0 * diameter),
        "minimum_norm_member": min_index,
        "diameter_pair": [int(first), int(second)],
        "float64_direct_vs_gram_signal_abs_diff": abs(signal - gram_signal),
        "float64_direct_vs_gram_diameter_abs_diff": abs(
            diameter - gram_diameter
        ),
        "classification_margin_ratio": (
            signal / diameter - 3.0 if diameter > 0 else float("inf")
        ),
    }


def evaluate(
    case_id: str,
    support_name: str,
    source_fidelity: str,
    cache: dict,
    dns: dict,
) -> tuple[list[dict], dict]:
    query_physical, query_transformed, weights, lx = J.prepare_support(
        case_id, source_fidelity, cache, dns
    )
    fields = {pathway: {} for pathway in J.PATHWAYS}
    for fidelity_id in J.FIDELITIES:
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
    for pathway in J.PATHWAYS:
        for fidelity_id in J.FIDELITIES:
            finite &= np.isfinite(fields[pathway][fidelity_id]).all(axis=1)
    finite &= np.isfinite(weights) & (weights > 0)
    weights = weights[finite]

    weighted_cells: dict[tuple[str, str], np.ndarray] = {}
    for pathway in J.PATHWAYS:
        for fidelity_id in J.FIDELITIES:
            field = fields[pathway][fidelity_id][finite]
            weighted_cells[(pathway, fidelity_id)] = weighted_vector(
                field, weights
            )

    CACHE_DIR.mkdir(exist_ok=True)
    vector_labels = [
        f"{pathway}:{fidelity_id}"
        for pathway in J.PATHWAYS
        for fidelity_id in J.FIDELITIES
    ]
    vector_matrix = np.stack(
        [
            weighted_cells[(pathway, fidelity_id)]
            for pathway in J.PATHWAYS
            for fidelity_id in J.FIDELITIES
        ]
    )
    cache_path = CACHE_DIR / f"{case_id}__{support_name}.npz"
    np.savez_compressed(
        cache_path,
        weighted_cell_vectors=vector_matrix,
        vector_labels=np.asarray(vector_labels),
        pathways=np.asarray(J.PATHWAYS),
        fidelities=np.asarray(J.FIDELITIES),
        coefficients=np.asarray(
            [J.COEFFICIENTS[fidelity_id] for fidelity_id in J.FIDELITIES],
            dtype=np.float64,
        ),
        finite_common_support_points=np.asarray([int(np.sum(finite))]),
        source_fidelity=np.asarray([source_fidelity]),
    )

    coherent_labels = list(J.PATHWAYS)
    coherent = []
    for pathway in coherent_labels:
        coherent.append(
            sum(
                J.COEFFICIENTS[fidelity_id]
                * weighted_cells[(pathway, fidelity_id)]
                for fidelity_id in J.FIDELITIES
            )
        )
    coherent_matrix = np.stack(coherent)

    recombination_labels = list(
        itertools.product(J.PATHWAYS, repeat=len(J.FIDELITIES))
    )
    recombined = []
    for choices in recombination_labels:
        recombined.append(
            sum(
                J.COEFFICIENTS[fidelity_id]
                * weighted_cells[(pathway, fidelity_id)]
                for fidelity_id, pathway in zip(J.FIDELITIES, choices)
            )
        )
    recombined_matrix = np.stack(recombined)

    coherent_metrics = set_metrics(coherent_matrix)
    recombined_metrics = set_metrics(recombined_matrix)
    rows = []
    for set_name, metrics in [
        ("coherent_four_pipeline_set", coherent_metrics),
        ("independent_4^4_recombination_set", recombined_metrics),
    ]:
        rows.append(
            {
                "case_id": case_id,
                "support": support_name,
                "source_fidelity": source_fidelity,
        "derived_vector_cache": str(cache_path),
                "finite_common_support_points": int(np.sum(finite)),
                "set_name": set_name,
                **metrics,
            }
        )

    details = {
        "case_id": case_id,
        "support": support_name,
        "source_fidelity": source_fidelity,
        "coherent_labels": coherent_labels,
        "recombination_fidelity_order": list(J.FIDELITIES),
        "recombination_pathway_order": [list(x) for x in recombination_labels],
        "coherent": coherent_metrics,
        "independent_recombination": recombined_metrics,
        "ratio_contraction": (
            coherent_metrics["ratio"] / recombined_metrics["ratio"]
            if recombined_metrics["ratio"] > 0
            else float("inf")
        ),
        "diameter_expansion": (
            recombined_metrics["diameter"] / coherent_metrics["diameter"]
            if coherent_metrics["diameter"] > 0
            else float("inf")
        ),
    }
    return rows, details


def main() -> None:
    cache = {
        (row["case_id"], row["fidelity_id"]): row
        for row in J.read_csv(J.INDEX)
    }
    dns = {
        (row["case_id"], row["fidelity_id"]): row
        for row in J.read_csv(J.DNS_INPUT)
    }
    rows = []
    details = []
    for case_id in J.PRIMARY_GEOMETRIES:
        for support_name, source_fidelity in (
            ("fine_F100", "F100"),
            ("shared_F000", "F000"),
        ):
            local_rows, local_details = evaluate(
                case_id, support_name, source_fidelity, cache, dns
            )
            rows.extend(local_rows)
            details.append(local_details)

    with (HERE / "07b_symmetric_mapping_recombination_audit.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    summary = {
        "protocol": "V12-20260730",
        "set_functionals": {
            "signal": "minimum weighted RMS norm over set members",
            "diameter": "maximum weighted RMS pairwise distance",
            "ratio": "signal / diameter",
            "threshold": "ratio >= 3 (frozen materiality stress rule)",
        },
        "results": details,
        "coherent_pass_count": sum(
            item["coherent"]["three_times_diameter_passed"]
            for item in details
        ),
        "independent_recombination_pass_count": sum(
            item["independent_recombination"][
                "three_times_diameter_passed"
            ]
            for item in details
        ),
        "interpretation": (
            "Both finite sets use identical signal, diameter, ratio, and "
            "threshold functionals. The comparison isolates cross-cell "
            "mapping coherence from the earlier triangle-envelope bound."
        ),
    }
    (HERE / "07b_symmetric_mapping_recombination_audit.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
