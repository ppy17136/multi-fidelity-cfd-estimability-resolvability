"""Run the self-contained V20 verification suite."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys
from pathlib import Path

sys.dont_write_bytecode = True

ROOT = Path(__file__).resolve().parent
TESTS = [
    "test_contrast_resolution.py",
    "test_mapping_set_metrics.py",
    "test_multilevel_closed_loop.py",
    "test_support_repair.py",
    "test_joint_mapping_ensemble.py",
    "test_second_physics_transfer.py",
]


def run(*args: str) -> None:
    subprocess.run(
        [sys.executable, *args],
        cwd=ROOT,
        check=True,
        env=dict(os.environ, PYTHONDONTWRITEBYTECODE="1"),
    )


def validate_archived_outputs() -> None:
    mapping = json.loads(
        (ROOT / "07b_symmetric_mapping_recombination_audit.json").read_text()
    )
    results = mapping["results"]
    assert len(results) == 4
    assert all(x["coherent"]["three_times_diameter_passed"] for x in results)
    assert all(
        not x["independent_recombination"]["three_times_diameter_passed"]
        for x in results
    )
    assert min(x["ratio_contraction"] for x in results) > 12.0
    loop = json.loads(
        (ROOT / "04b_multilevel_closed_loop_summary.json").read_text()
    )
    assert loop["protocol"] == "V20-20260803"
    assert loop["paired_action_indexed_scenarios"] is True
    assert len(loop["policies"]) == 7
    assert len(loop["policies"]) == len(set(loop["policies"]))
    assert "c_optimal" not in loop["policies"]
    grid = json.loads(
        (ROOT / "04c_truth_grid_closed_loop_summary.json").read_text()
    )
    assert grid["protocol"] == "V20-20260803"
    assert grid["paired_action_indexed_scenarios"] is True
    assert grid["policies"] == loop["policies"]
    assert grid["all_wrong_decisive_fractions_zero"] is True
    assert len(grid["true_deltas"]) == 11


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--rebuild-synthetic", action="store_true")
    args = parser.parse_args()
    run("-m", "pytest", "-q", "-p", "no:cacheprovider", *TESTS)
    validate_archived_outputs()
    if args.rebuild_synthetic:
        for script in (
            "01_run_support_repair_benchmarks.py",
            "01b_run_randomized_support_benchmarks.py",
            "02_run_gate2_calibration.py",
            "04b_run_multilevel_closed_loop_acquisition.py",
            "04c_run_truth_grid_closed_loop_audit.py",
        ):
            run(script)
    run("07c_recompute_mapping_from_derived_vectors.py", "--check-only")
    run(
        "V16_complete_lattice_mapping_audit/run_v16_audit.py",
        "--check-only",
    )
    print("V20 verification passed.")


if __name__ == "__main__":
    main()
