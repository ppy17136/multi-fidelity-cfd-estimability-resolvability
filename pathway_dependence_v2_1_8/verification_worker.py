"""One-command verification for the pathway-coupling certificate release."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def run(*args: str) -> None:
    print("+", " ".join(args), flush=True)
    subprocess.run(args, cwd=ROOT, check=True)


def require_close(actual: float, expected: float, tol: float, label: str) -> None:
    if abs(actual - expected) > tol:
        raise RuntimeError(f"{label}: observed {actual!r}, expected {expected!r}")


def main() -> None:
    run(sys.executable, "-m", "unittest", "discover", "-v")
    run(sys.executable, "run_real_cfd_common_partition_profiles.py")
    run(sys.executable, "build_authoritative_cfd_summary.py")
    if "--include-robust-union" in sys.argv:
        print("AUXILIARY ROBUST-UNION RECONSTRUCTION (not main CFD claims)",flush=True)
        for name in ("run_robust_union_profiles.py", "analyze_robust_union_epsilon.py",
                     "build_robust_union_threshold_frontiers.py", "run_robust_union_edge_cost_sensitivity.py"):
            run(sys.executable, "auxiliary/robust_union/"+name)

    summary = json.loads((ROOT / "42_REAL_CFD_COMMON_PARTITION_SUMMARY.json").read_text(encoding="utf-8"))
    if summary["audit_count"] != 8:
        raise RuntimeError("Expected eight CFD mapping audits")
    interval = summary["all_audit_common_threshold_interval"]
    if not interval["nonempty"]:
        raise RuntimeError("Expected a nonempty common threshold interval")
    require_close(interval["lower_open"], 1.6442596248017918, 1e-12, "interval lower endpoint")
    require_close(interval["upper_closed"], 2.2831676147449826, 1e-12, "interval upper endpoint")
    epsilon_costs = sorted({row["D_cost_at_relative_tolerance_1e-7"] for row in summary["rows"]})
    if epsilon_costs != [0.5]:
        raise RuntimeError("All balanced-cost near-diameter certificates must equal 0.5")
    exact_cost_counts = {
        "0.5": sum(row["exact_independent_D_cost"] == 0.5 for row in summary["rows"]),
        "1.0": sum(row["exact_independent_D_cost"] == 1.0 for row in summary["rows"]),
    }
    if exact_cost_counts != {"0.5": 4, "1.0": 4}:
        raise RuntimeError("Unexpected exact floating-array diameter cost split")

    run(sys.executable, "verify_worked_example.py")
    run(sys.executable, "analyze_weight_stability.py")
    run(sys.executable, "verify_algorithm_table.py")
    run(sys.executable, "verification_extensions/test_new_theory.py")
    run(sys.executable, "stochastic_diffusion/verify_certificates.py")
    run(sys.executable, "stochastic_diffusion/verify_reference_solver.py")
    run(sys.executable, "verification_extensions/test_readonly_verification.py")

    run(sys.executable, "verification_extensions/test_release_semantics.py")
    print("Isolated reconstruction and witness checks passed")


if __name__ == "__main__":
    main()