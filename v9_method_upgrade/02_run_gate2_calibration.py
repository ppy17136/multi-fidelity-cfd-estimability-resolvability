"""Frozen Gate-2 implementation verification and decision calibration."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from contrast_resolution import (
    contrast_standard_deviation,
    deterministic_contrast_bound,
    deterministic_resolution_decision,
    legacy_ratio_decision,
    probabilistic_resolution_decision,
)


ROOT = Path(__file__).resolve().parent
SEED = 20260730
COEFFICIENTS = np.asarray([1.0, -1.0, -1.0, 1.0])
TRUE_CELL_BOUND = 0.25
TRUE_CONTRAST_BOUND = deterministic_contrast_bound(
    COEFFICIENTS, np.full(4, TRUE_CELL_BOUND)
)
MINIMUM_EFFECT = 3.0
N = 200_000


def independent_bounded_errors(rng: np.random.Generator, n: int) -> np.ndarray:
    return rng.uniform(-TRUE_CELL_BOUND, TRUE_CELL_BOUND, size=(n, 4))


def correlated_normal_errors(
    rng: np.random.Generator, n: int, correlation: float
) -> tuple[np.ndarray, np.ndarray]:
    sd = TRUE_CELL_BOUND
    covariance = np.full((4, 4), correlation * sd * sd)
    np.fill_diagonal(covariance, sd * sd)
    return rng.multivariate_normal(np.zeros(4), covariance, size=n), covariance


def classify(true_delta: float) -> bool:
    return abs(true_delta) > MINIMUM_EFFECT


def summarize_binary_decisions(
    true_delta: float, resolved: np.ndarray
) -> dict[str, float | int | bool]:
    truth = classify(true_delta)
    resolved = np.asarray(resolved, dtype=bool)
    false_resolved = np.mean(resolved & (not truth))
    false_unresolved = np.mean((~resolved) & truth)
    return {
        "true_delta": true_delta,
        "scientifically_above_minimum_effect": truth,
        "resolved_fraction": float(np.mean(resolved)),
        "false_resolved_fraction": float(false_resolved),
        "false_unresolved_fraction": float(false_unresolved),
        "correct_decision_fraction": float(
            np.mean(resolved == np.full(len(resolved), truth))
        ),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)
    bounded = independent_bounded_errors(rng, N)
    contrast_error = bounded @ COEFFICIENTS
    deltas = [0.0, 2.0, 3.0, 5.0]
    scenario_rows: list[dict] = []

    for delta in deltas:
        observed = delta + contrast_error
        for threshold in np.arange(1.5, 5.01, 0.25):
            resolved = np.abs(observed) / TRUE_CONTRAST_BOUND > threshold
            row = summarize_binary_decisions(delta, resolved)
            row |= {
                "error_model": "independent_bounded",
                "decision_method": "legacy_ratio",
                "threshold_or_confidence": float(threshold),
                "assumed_floor_fraction": 1.0,
            }
            scenario_rows.append(row)

        deterministic = np.asarray(
            [
                deterministic_resolution_decision(
                    value, TRUE_CONTRAST_BOUND, MINIMUM_EFFECT
                ).resolved
                for value in observed
            ]
        )
        row = summarize_binary_decisions(delta, deterministic)
        row |= {
            "error_model": "independent_bounded",
            "decision_method": "deterministic_lower_bound",
            "threshold_or_confidence": None,
            "assumed_floor_fraction": 1.0,
        }
        scenario_rows.append(row)

        for assumed_fraction in (0.9, 0.75, 0.5):
            assumed_bound = TRUE_CONTRAST_BOUND * assumed_fraction
            underestimated = np.asarray(
                [
                    deterministic_resolution_decision(
                        value, assumed_bound, MINIMUM_EFFECT
                    ).resolved
                    for value in observed
                ]
            )
            row = summarize_binary_decisions(delta, underestimated)
            row |= {
                "error_model": "independent_bounded",
                "decision_method": "deterministic_lower_bound",
                "threshold_or_confidence": None,
                "assumed_floor_fraction": assumed_fraction,
            }
            scenario_rows.append(row)

    covariance_rows: list[dict] = []
    for correlation in (-0.25, 0.0, 0.5, 0.9):
        errors, covariance = correlated_normal_errors(rng, N, correlation)
        contrast_error_normal = errors @ COEFFICIENTS
        contrast_sd = contrast_standard_deviation(COEFFICIENTS, covariance)
        for delta in deltas:
            observed = delta + contrast_error_normal
            resolved = np.asarray(
                [
                    probabilistic_resolution_decision(
                        value,
                        contrast_sd,
                        MINIMUM_EFFECT,
                        confidence=0.95,
                    ).resolved
                    for value in observed
                ]
            )
            row = summarize_binary_decisions(delta, resolved)
            row |= {
                "error_model": "correlated_normal",
                "correlation": correlation,
                "contrast_sd": contrast_sd,
                "decision_method": "normal_covariance",
                "threshold_or_confidence": 0.95,
                "assumed_floor_fraction": 1.0,
            }
            covariance_rows.append(row)

    all_rows = scenario_rows + covariance_rows
    fieldnames = sorted({key for row in all_rows for key in row})
    with (ROOT / "02_gate2_calibration_results.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(all_rows)

    primary = [
        row
        for row in scenario_rows
        if row["decision_method"] == "deterministic_lower_bound"
        and row["assumed_floor_fraction"] == 1.0
    ]
    normal_primary = [
        row
        for row in covariance_rows
        if row["threshold_or_confidence"] == 0.95
    ]
    summary = {
        "protocol": "V9-20260730",
        "seed": SEED,
        "realisations_per_scenario": N,
        "true_contrast_bound": TRUE_CONTRAST_BOUND,
        "minimum_effect": MINIMUM_EFFECT,
        "deterministic_primary": primary,
        "normal_covariance_primary": normal_primary,
        "maximum_false_resolved_normal_95": max(
            row["false_resolved_fraction"] for row in normal_primary
        ),
        "maximum_false_unresolved_normal_95": max(
            row["false_unresolved_fraction"] for row in normal_primary
        ),
        "claim_boundary": (
            "Normal-model error rates apply only to the constructed covariance "
            "experiment; deterministic CFD sensitivities remain nonprobabilistic."
        ),
    }
    (ROOT / "02_gate2_calibration_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

