"""Recompute the symmetric mapping audit from released float64 vector caches."""

from __future__ import annotations

import argparse
import itertools
import json
from pathlib import Path

import numpy as np

from mapping_set_metrics import set_metrics


HERE = Path(__file__).resolve().parent
CACHE_DIR = HERE / "derived_mapping_vectors_v12"
REFERENCE = HERE / "07b_symmetric_mapping_recombination_audit.json"
OUTPUT = HERE / "07c_recomputed_mapping_from_vectors.json"


def evaluate(path: Path) -> dict:
    with np.load(path) as data:
        vectors = np.asarray(data["weighted_cell_vectors"], dtype=np.float64)
        pathways = [str(x) for x in data["pathways"]]
        fidelities = [str(x) for x in data["fidelities"]]
        coefficients = np.asarray(data["coefficients"], dtype=np.float64)
        point_count = int(data["finite_common_support_points"][0])

    expected_shape = (len(pathways) * len(fidelities), 2 * point_count)
    if vectors.shape != expected_shape:
        raise RuntimeError(
            f"{path.name}: vector shape {vectors.shape} != {expected_shape}"
        )
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
    return {
        "cache": path.name,
        "finite_common_support_points": point_count,
        "coherent": set_metrics(coherent),
        "independent_recombination": set_metrics(recombined),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--check-only", action="store_true")
    args = parser.parse_args()
    results = [evaluate(path) for path in sorted(CACHE_DIR.glob("*.npz"))]
    if len(results) != 4:
        raise RuntimeError(f"expected four vector caches, found {len(results)}")

    reference = json.loads(REFERENCE.read_text(encoding="utf-8"))
    reference_by_key = {
        f"{item['case_id']}__{item['support']}.npz": item
        for item in reference["results"]
    }
    maximum_abs_difference = 0.0
    for item in results:
        ref = reference_by_key[item["cache"]]
        for set_name in ("coherent", "independent_recombination"):
            for metric in ("minimum_norm", "diameter", "ratio"):
                difference = abs(
                    float(item[set_name][metric])
                    - float(ref[set_name][metric])
                )
                maximum_abs_difference = max(
                    maximum_abs_difference, difference
                )
                if difference > 1e-12:
                    raise RuntimeError(
                        f"{item['cache']} {set_name} {metric}: "
                        f"difference {difference:.3e}"
                    )

    output = {
        "protocol": "V12-20260730",
        "input": "released float64 weighted cell-vector caches",
        "result_count": len(results),
        "maximum_abs_difference_from_primary_audit": maximum_abs_difference,
        "passed": True,
        "results": results,
    }
    if not args.check_only:
        OUTPUT.write_text(
            json.dumps(output, indent=2) + "\n", encoding="utf-8"
        )
    print(json.dumps(output, indent=2))


if __name__ == "__main__":
    main()
