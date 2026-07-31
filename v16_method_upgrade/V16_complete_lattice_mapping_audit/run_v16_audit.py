"""Reconstruct and verify all eight V16 mapping-dependence audits.

The script uses only files under this directory. No author workstation path,
OpenFOAM case, historical version directory, or unpublished cache is required.
"""

from __future__ import annotations

import argparse
import csv
import itertools
import json
from pathlib import Path

import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "input"
EXPECTED = HERE / "expected"
GENERATED = HERE / "generated"
CASES = (
    "alph05-10071-2024",
    "alph15-10929-3036",
    "alph15-13929-2024",
    "alph15-7929-2024",
)
SUPPORTS = ("fine_F100", "shared_F000")
THRESHOLD = 3.0
TOLERANCE = 1.0e-12


def set_metrics(vectors: np.ndarray) -> dict:
    x = np.asarray(vectors, dtype=np.float64)
    gram = x @ x.T
    norm2 = np.maximum(np.diag(gram), 0.0)
    distance2 = norm2[:, None] + norm2[None, :] - 2.0 * gram
    np.maximum(distance2, 0.0, out=distance2)
    min_index = int(np.argmin(norm2))
    maximum_flat = int(np.argmax(distance2))
    first, second = np.unravel_index(maximum_flat, distance2.shape)
    signal = float(np.linalg.norm(x[min_index]))
    diameter = float(np.linalg.norm(x[first] - x[second]))
    return {
        "set_size": int(len(x)),
        "minimum_norm": signal,
        "diameter": diameter,
        "ratio": signal / diameter if diameter > 0 else float("inf"),
        "three_times_diameter_passed": bool(signal >= THRESHOLD * diameter),
        "minimum_norm_member": min_index,
        "diameter_pair": [int(first), int(second)],
    }


def evaluate(path: Path) -> tuple[list[dict], dict]:
    with np.load(path) as data:
        vectors = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(x) for x in data["pathways"]]
        fidelities = [str(x) for x in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        point_count = int(data["finite_common_support_points"][0])
        case_id = str(data["case_id"][0])
        support = str(data["support"][0])
        source_fidelity = str(data["source_fidelity"][0])

    if path.name != f"{case_id}__{support}.npz":
        raise RuntimeError(f"{path.name}: embedded case/support do not match name")
    expected_shape = (len(pathways) * len(fidelities), 2 * point_count)
    if vectors.shape != expected_shape or vectors.dtype != np.float64:
        raise RuntimeError(
            f"{path.name}: vectors {vectors.shape}/{vectors.dtype}, "
            f"expected {expected_shape}/float64"
        )
    if len(pathways) != 4 or len(fidelities) != 4:
        raise RuntimeError(f"{path.name}: expected a four-by-four pathway lattice")

    cell = {
        (pathway, fidelity): vectors[p * len(fidelities) + f]
        for p, pathway in enumerate(pathways)
        for f, fidelity in enumerate(fidelities)
    }
    coherent = np.stack(
        [
            sum(
                coefficients[f] * cell[(pathway, fidelity)]
                for f, fidelity in enumerate(fidelities)
            )
            for pathway in pathways
        ]
    )
    combinations = list(itertools.product(pathways, repeat=len(fidelities)))
    recombined = np.stack(
        [
            sum(
                coefficients[f] * cell[(pathway, fidelity)]
                for f, (fidelity, pathway) in enumerate(
                    zip(fidelities, choices)
                )
            )
            for choices in combinations
        ]
    )
    coherent_metrics = set_metrics(coherent)
    recombined_metrics = set_metrics(recombined)
    rows = []
    for set_name, metrics in (
        ("coherent_four_pipeline_set", coherent_metrics),
        ("independent_4^4_recombination_set", recombined_metrics),
    ):
        rows.append(
            {
                "case_id": case_id,
                "support": support,
                "source_fidelity": source_fidelity,
                "finite_common_support_points": point_count,
                "set_name": set_name,
                **metrics,
            }
        )
    detail = {
        "case_id": case_id,
        "support": support,
        "source_fidelity": source_fidelity,
        "coherent_labels": pathways,
        "recombination_fidelity_order": fidelities,
        "coherent": coherent_metrics,
        "independent_recombination": recombined_metrics,
        "ratio_contraction": (
            coherent_metrics["ratio"] / recombined_metrics["ratio"]
        ),
        "diameter_expansion": (
            recombined_metrics["diameter"] / coherent_metrics["diameter"]
        ),
    }
    return rows, detail


def check_eligibility_index() -> None:
    with (INPUT / "complete_lattice_index.csv").open(
        newline="", encoding="utf-8"
    ) as handle:
        rows = list(csv.DictReader(handle))
    eligible = tuple(
        row["case_id"]
        for row in rows
        if row["eligible_complete_lattice"].lower() == "true"
    )
    included = tuple(
        row["case_id"]
        for row in rows
        if row["included_in_eight_audit"].lower() == "true"
    )
    if eligible != CASES or included != CASES:
        raise RuntimeError(
            f"eligibility/inclusion mismatch: eligible={eligible}, "
            f"included={included}, expected={CASES}"
        )


def compare_metric(name: str, actual: float, expected: float) -> float:
    difference = abs(float(actual) - float(expected))
    if difference > TOLERANCE:
        raise RuntimeError(
            f"{name}: absolute difference {difference:.3e} exceeds {TOLERANCE:.1e}"
        )
    return difference


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--check-only",
        action="store_true",
        help="verify without writing generated outputs",
    )
    args = parser.parse_args()

    check_eligibility_index()
    expected_names = {
        f"{case_id}__{support}.npz"
        for case_id in CASES
        for support in SUPPORTS
    }
    observed_names = {path.name for path in INPUT.glob("*.npz")}
    if observed_names != expected_names:
        raise RuntimeError(
            f"vector cache names differ: {sorted(observed_names ^ expected_names)}"
        )

    rows: list[dict] = []
    details: list[dict] = []
    for case_id in CASES:
        for support in SUPPORTS:
            local_rows, detail = evaluate(INPUT / f"{case_id}__{support}.npz")
            rows.extend(local_rows)
            details.append(detail)

    coherent = [item["coherent"] for item in details]
    recombined = [item["independent_recombination"] for item in details]
    coherent_min = min(item["ratio"] for item in coherent)
    recombined_max = max(item["ratio"] for item in recombined)
    summary = {
        "protocol": "V16-ALL-COMPLETE-LATTICES-20260731",
        "objective_inclusion_rule_cases": list(CASES),
        "supports_per_case": list(SUPPORTS),
        "audit_count": len(details),
        "coherent_set_size": 4,
        "independent_set_size": 256,
        "numeric_precision": "float64 weighted vectors and Gram matrices",
        "frozen_display_threshold": THRESHOLD,
        "coherent_pass_count": sum(
            item["three_times_diameter_passed"] for item in coherent
        ),
        "independent_pass_count": sum(
            item["three_times_diameter_passed"] for item in recombined
        ),
        "paired_classification_reversal_count": sum(
            c["three_times_diameter_passed"]
            and not r["three_times_diameter_passed"]
            for c, r in zip(coherent, recombined)
        ),
        "directional_ratio_reduction_count": sum(
            r["ratio"] < c["ratio"] for c, r in zip(coherent, recombined)
        ),
        "coherent_ratio_minimum": coherent_min,
        "independent_ratio_maximum": recombined_max,
        "common_separation_threshold_interval": [
            recombined_max,
            coherent_min,
        ],
        "diameter_expansion_range": [
            min(item["diameter_expansion"] for item in details),
            max(item["diameter_expansion"] for item in details),
        ],
        "ratio_contraction_range": [
            min(item["ratio_contraction"] for item in details),
            max(item["ratio_contraction"] for item in details),
        ],
        "results": details,
        "interpretation": (
            "Every lattice meeting the frozen completeness rule was evaluated. "
            "This is a mapping-dependence classification, not a total-CFD "
            "uncertainty or complete Gate-2 decision."
        ),
    }

    if summary["audit_count"] != 8:
        raise RuntimeError("expected eight audits")
    if summary["coherent_pass_count"] != 7:
        raise RuntimeError("expected seven coherent passes at R=3")
    if summary["independent_pass_count"] != 0:
        raise RuntimeError("expected zero independent-recombination passes at R=3")
    if summary["paired_classification_reversal_count"] != 7:
        raise RuntimeError("expected seven paired reversals at R=3")
    if summary["directional_ratio_reduction_count"] != 8:
        raise RuntimeError("expected ratio reduction in all eight audit pairs")
    if not recombined_max < coherent_min:
        raise RuntimeError("common separation interval is empty")

    expected = json.loads((EXPECTED / "audit.json").read_text(encoding="utf-8"))
    expected_by_key = {
        (item["case_id"], item["support"]): item for item in expected["results"]
    }
    maximum_abs_difference = 0.0
    for item in details:
        key = (item["case_id"], item["support"])
        reference = expected_by_key[key]
        for set_name in ("coherent", "independent_recombination"):
            for metric in ("minimum_norm", "diameter", "ratio"):
                maximum_abs_difference = max(
                    maximum_abs_difference,
                    compare_metric(
                        f"{key}/{set_name}/{metric}",
                        item[set_name][metric],
                        reference[set_name][metric],
                    ),
                )
        for metric in ("ratio_contraction", "diameter_expansion"):
            maximum_abs_difference = max(
                maximum_abs_difference,
                compare_metric(
                    f"{key}/{metric}",
                    item[metric],
                    reference[metric],
                ),
            )
    for metric in (
        "coherent_ratio_minimum",
        "independent_ratio_maximum",
    ):
        maximum_abs_difference = max(
            maximum_abs_difference,
            compare_metric(metric, summary[metric], expected[metric]),
        )

    verification = {
        **summary,
        "maximum_abs_difference_from_expected": maximum_abs_difference,
        "passed": True,
    }
    if not args.check_only:
        GENERATED.mkdir(exist_ok=True)
        write_csv(GENERATED / "audit.csv", rows)
        (GENERATED / "audit.json").write_text(
            json.dumps(verification, indent=2) + "\n", encoding="utf-8"
        )
        (GENERATED / "V16_EIGHT_AUDIT_RECONSTRUCTION.ok").write_text(
            "V16 portable eight-audit reconstruction passed.\n",
            encoding="utf-8",
        )
    print(
        json.dumps(
            {
                "protocol": verification["protocol"],
                "audit_count": verification["audit_count"],
                "coherent_pass_count": verification["coherent_pass_count"],
                "independent_pass_count": verification["independent_pass_count"],
                "paired_classification_reversal_count": verification[
                    "paired_classification_reversal_count"
                ],
                "directional_ratio_reduction_count": verification[
                    "directional_ratio_reduction_count"
                ],
                "common_separation_threshold_interval": verification[
                    "common_separation_threshold_interval"
                ],
                "maximum_abs_difference_from_expected": verification[
                    "maximum_abs_difference_from_expected"
                ],
                "passed": verification["passed"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
