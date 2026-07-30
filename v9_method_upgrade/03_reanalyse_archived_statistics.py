"""Retrospective uncertainty intervals and threshold sensitivity.

All analyses in this file use already opened outcomes and are explicitly
retrospective.  Bootstrap resampling is performed at the independent replicate
level rather than treating learner pairs as independent experimental units.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parent
PROJECT = ROOT.parent
SEED = 20260730
BOOTSTRAPS = 100_000


def percentile_interval(values: np.ndarray) -> dict[str, float]:
    q = np.quantile(values, [0.025, 0.5, 0.975])
    return {"lower_95": float(q[0]), "median": float(q[1]), "upper_95": float(q[2])}


def bootstrap_median(
    values: np.ndarray, indices: np.ndarray
) -> dict[str, float]:
    samples = np.median(values[indices], axis=1)
    return percentile_interval(samples)


def wilson_interval(successes: int, n: int, z: float = 1.959963984540054) -> dict:
    p = successes / n
    denominator = 1.0 + z * z / n
    centre = (p + z * z / (2.0 * n)) / denominator
    half = (
        z
        * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n))
        / denominator
    )
    return {
        "successes": successes,
        "trials": n,
        "fraction": p,
        "lower_95": float(centre - half),
        "upper_95": float(centre + half),
    }


def rashomon_intervals(rng: np.random.Generator) -> dict:
    data = pd.read_csv(PROJECT / "166_rashomon_replicate_summary.csv")
    n = len(data)
    indices = rng.integers(0, n, size=(BOOTSTRAPS, n))
    aggregate = data["median_aggregate_prediction_disagreement_nrmse"].to_numpy()
    interaction = data[
        "median_interaction_prediction_disagreement_nrmse"
    ].to_numpy()
    within_ratio = data[
        "median_interaction_to_aggregate_disagreement_ratio"
    ].to_numpy()
    ratio_of_medians_samples = np.median(interaction[indices], axis=1) / np.median(
        aggregate[indices], axis=1
    )
    return {
        "independent_replicates": n,
        "learner_pairs_total": 480,
        "aggregate_disagreement_median": float(np.median(aggregate)),
        "aggregate_disagreement_bootstrap_95": bootstrap_median(aggregate, indices),
        "interaction_disagreement_median": float(np.median(interaction)),
        "interaction_disagreement_bootstrap_95": bootstrap_median(
            interaction, indices
        ),
        "median_within_pair_ratio": float(np.median(within_ratio)),
        "median_within_pair_ratio_bootstrap_95": bootstrap_median(
            within_ratio, indices
        ),
        "ratio_of_reported_medians": float(
            np.median(interaction) / np.median(aggregate)
        ),
        "ratio_of_medians_bootstrap_95": percentile_interval(
            ratio_of_medians_samples
        ),
        "wording": (
            "The interaction and aggregate medians and the median within-pair "
            "ratio are distinct statistics."
        ),
    }


def equal_cost_intervals(rng: np.random.Generator) -> dict:
    data = pd.read_csv(PROJECT / "168_equal_cost_acquisition_unit_metrics.csv")
    pivot = data.pivot(
        index="replicate", columns="policy", values="median_interaction_nrmse"
    )
    baseline = pivot["F110_heavy"].to_numpy()
    proposed = pivot["contrast_complete"].to_numpy()
    n = len(pivot)
    indices = rng.integers(0, n, size=(BOOTSTRAPS, n))
    ratio_of_medians_reduction = 1.0 - np.median(
        proposed[indices], axis=1
    ) / np.median(baseline[indices], axis=1)
    paired_relative = (baseline - proposed) / baseline
    paired_difference = baseline - proposed
    wins = int(np.sum(proposed < baseline))
    return {
        "independent_replicates": n,
        "original_ratio_of_medians_reduction": float(
            1.0 - np.median(proposed) / np.median(baseline)
        ),
        "ratio_of_medians_reduction_bootstrap_95": percentile_interval(
            ratio_of_medians_reduction
        ),
        "paired_relative_reduction_median": float(np.median(paired_relative)),
        "paired_relative_reduction_bootstrap_95": bootstrap_median(
            paired_relative, indices
        ),
        "paired_absolute_difference_median": float(np.median(paired_difference)),
        "paired_absolute_difference_bootstrap_95": bootstrap_median(
            paired_difference, indices
        ),
        "win_probability_wilson_95": wilson_interval(wins, n),
    }


def cfd_threshold_sensitivity() -> tuple[list[dict], dict]:
    source = json.loads(
        (PROJECT / "188_exploratory_four_lattice_materiality.json").read_text(
            encoding="utf-8"
        )
    )
    ratios = {
        row["case_id"]: float(row["exact_to_archived_floor_ratio"])
        for row in source["per_geometry"]
    }
    ratios["pooled_four"] = float(
        source["pooled_all_four_archived_floor"]["exact_to_mapping_floor_ratio"]
    )
    ratios["pooled_four_conservative"] = float(
        source["pooled_all_four_conservative_floor"]["exact_to_mapping_floor_ratio"]
    )
    rows = []
    for threshold in np.arange(1.5, 5.01, 0.25):
        row = {"threshold": float(threshold)}
        for name, ratio in ratios.items():
            row[name] = bool(ratio > threshold)
        row["individual_cases_passing"] = int(
            sum(ratios[name] > threshold for name in ratios if name.startswith("alph"))
        )
        rows.append(row)
    summary = {
        "ratios": ratios,
        "thresholds_at_which_individual_decisions_change": sorted(
            float(value) for key, value in ratios.items() if key.startswith("alph")
        ),
        "frozen_primary_threshold": 3.0,
        "all_individual_unresolved_at_frozen_threshold": all(
            value <= 3.0 for key, value in ratios.items() if key.startswith("alph")
        ),
        "claim_boundary": (
            "This retrospective threshold curve does not replace the frozen "
            "primary decision."
        ),
    }
    return rows, summary


def main() -> None:
    rng = np.random.default_rng(SEED)
    rashomon = rashomon_intervals(rng)
    equal_cost = equal_cost_intervals(rng)
    sensitivity_rows, cfd_summary = cfd_threshold_sensitivity()

    with (ROOT / "03_cfd_threshold_sensitivity.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=list(sensitivity_rows[0]))
        writer.writeheader()
        writer.writerows(sensitivity_rows)

    summary = {
        "status": "retrospective_reanalysis",
        "seed": SEED,
        "bootstrap_replicates": BOOTSTRAPS,
        "rashomon": rashomon,
        "equal_cost_acquisition": equal_cost,
        "cfd_threshold_sensitivity": cfd_summary,
    }
    (ROOT / "03_retrospective_statistics_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()

