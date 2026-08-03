"""Frozen Gate-2 three-state implementation verification and calibration."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

from contrast_resolution import (
    contrast_standard_deviation,
    deterministic_contrast_bound,
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
DELTAS = [0.0, 1.0, 2.0, 2.5, 2.9, 3.0, 3.1, 3.5, 4.0, 5.0, 6.0]


def independent_bounded_errors(rng: np.random.Generator, n: int) -> np.ndarray:
    return rng.uniform(-TRUE_CELL_BOUND, TRUE_CELL_BOUND, size=(n, 4))


def correlated_normal_errors(
    rng: np.random.Generator, n: int, correlation: float
) -> tuple[np.ndarray, np.ndarray]:
    sd = TRUE_CELL_BOUND
    covariance = np.full((4, 4), correlation * sd * sd)
    np.fill_diagonal(covariance, sd * sd)
    return rng.multivariate_normal(np.zeros(4), covariance, size=n), covariance


def truth_state(true_delta: float) -> str:
    magnitude = abs(true_delta)
    if magnitude > MINIMUM_EFFECT:
        return "effect_present"
    if magnitude < MINIMUM_EFFECT:
        return "below_minimum_effect"
    return "boundary"


def interval_states(
    estimates: np.ndarray, radius: float, minimum_effect: float
) -> np.ndarray:
    magnitude = np.abs(estimates)
    lower = np.maximum(magnitude - radius, 0.0)
    upper = magnitude + radius
    states = np.full(len(estimates), "indeterminate", dtype=object)
    states[lower > minimum_effect] = "effect_present"
    states[upper < minimum_effect] = "below_minimum_effect"
    return states


def summarize_decision_states(
    true_delta: float, states: np.ndarray
) -> dict[str, float | str]:
    truth = truth_state(true_delta)
    states = np.asarray(states, dtype=object)
    present = states == "effect_present"
    below = states == "below_minimum_effect"
    indeterminate = states == "indeterminate"
    if truth == "effect_present":
        wrong_decisive = below
        correct_decisive = present
    elif truth == "below_minimum_effect":
        wrong_decisive = present
        correct_decisive = below
    else:
        wrong_decisive = present | below
        correct_decisive = np.zeros(len(states), dtype=bool)
    return {
        "true_delta": true_delta,
        "truth_state": truth,
        "effect_present_fraction": float(np.mean(present)),
        "below_minimum_effect_fraction": float(np.mean(below)),
        "indeterminate_fraction": float(np.mean(indeterminate)),
        "wrong_decisive_fraction": float(np.mean(wrong_decisive)),
        "correct_decisive_fraction": float(np.mean(correct_decisive)),
    }


def main() -> None:
    rng = np.random.default_rng(SEED)
    bounded = independent_bounded_errors(rng, N)
    contrast_error = bounded @ COEFFICIENTS
    scenario_rows: list[dict] = []

    for delta in DELTAS:
        observed = delta + contrast_error
        for threshold in np.arange(1.5, 5.01, 0.25):
            states = np.where(
                np.abs(observed) / TRUE_CONTRAST_BOUND > threshold,
                "effect_present",
                "indeterminate",
            )
            row = summarize_decision_states(delta, states)
            row |= {
                "error_model": "independent_bounded",
                "decision_method": "legacy_ratio",
                "threshold_or_confidence": float(threshold),
                "assumed_floor_fraction": 1.0,
            }
            scenario_rows.append(row)

        states = interval_states(
            observed, TRUE_CONTRAST_BOUND, MINIMUM_EFFECT
        )
        row = summarize_decision_states(delta, states)
        row |= {
            "error_model": "independent_bounded",
            "decision_method": "deterministic_interval",
            "threshold_or_confidence": None,
            "assumed_floor_fraction": 1.0,
        }
        scenario_rows.append(row)

        for assumed_fraction in (0.9, 0.75, 0.5):
            states = interval_states(
                observed,
                TRUE_CONTRAST_BOUND * assumed_fraction,
                MINIMUM_EFFECT,
            )
            row = summarize_decision_states(delta, states)
            row |= {
                "error_model": "independent_bounded",
                "decision_method": "deterministic_interval",
                "threshold_or_confidence": None,
                "assumed_floor_fraction": assumed_fraction,
            }
            scenario_rows.append(row)

    covariance_rows: list[dict] = []
    for correlation in (-0.25, 0.0, 0.5, 0.9):
        errors, covariance = correlated_normal_errors(rng, N, correlation)
        contrast_error_normal = errors @ COEFFICIENTS
        contrast_sd = contrast_standard_deviation(COEFFICIENTS, covariance)
        radius = 1.959963984540054 * contrast_sd
        for delta in DELTAS:
            observed = delta + contrast_error_normal
            states = interval_states(observed, radius, MINIMUM_EFFECT)
            row = summarize_decision_states(delta, states)
            row |= {
                "error_model": "correlated_normal",
                "correlation": correlation,
                "contrast_sd": contrast_sd,
                "decision_method": "central_normal_interval",
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
        if row["decision_method"] == "deterministic_interval"
        and row["assumed_floor_fraction"] == 1.0
    ]
    normal_primary = [
        row
        for row in covariance_rows
        if row["threshold_or_confidence"] == 0.95
    ]
    summary = {
        "protocol": "V12-20260730", "frozen_parent_protocol": "V10-20260730",
        "seed": SEED,
        "realisations_per_scenario": N,
        "tested_true_contrasts": DELTAS,
        "true_contrast_bound": TRUE_CONTRAST_BOUND,
        "minimum_effect": MINIMUM_EFFECT,
        "decision_semantics": (
            "effect_present if the lower interval endpoint exceeds the "
            "minimum effect; below_minimum_effect if the upper endpoint is "
            "below it; otherwise indeterminate"
        ),
        "deterministic_primary": primary,
        "normal_covariance_primary": normal_primary,
        "maximum_wrong_decisive_normal_95": max(
            row["wrong_decisive_fraction"] for row in normal_primary
        ),
        "normal_interval": "central two-sided 95% normal interval",
        "claim_boundary": (
            "Normal-model error rates apply only to the constructed covariance "
            "experiment; deterministic CFD sensitivities remain "
            "nonprobabilistic."
        ),
    }
    (ROOT / "02_gate2_calibration_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
