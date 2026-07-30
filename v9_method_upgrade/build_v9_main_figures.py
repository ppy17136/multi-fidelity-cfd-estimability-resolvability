from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch


HERE = Path(__file__).resolve().parent
ROOT = HERE
OUT = HERE.parent / "figures" / "v9_method_upgrade"
OUT.mkdir(parents=True, exist_ok=True)

plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 9.5,
        "axes.titlesize": 11,
        "axes.labelsize": 9.5,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    }
)

NAVY = "#1F3A5F"
BLUE = "#2F6B9A"
CYAN = "#67A9CF"
ORANGE = "#E07A3F"
RED = "#B4473A"
GREEN = "#3D8B74"
PURPLE = "#7353A4"
GREY = "#6B7280"
LIGHT = "#F4F7FA"


def save(fig: plt.Figure, stem: str) -> None:
    fig.savefig(OUT / f"{stem}.png", bbox_inches="tight", facecolor="white")
    fig.savefig(OUT / f"{stem}.pdf", bbox_inches="tight", facecolor="white")
    plt.close(fig)


def box(ax, xy, width, height, text, edge, face="white", fontsize=9.2):
    patch = FancyBboxPatch(
        xy,
        width,
        height,
        boxstyle="round,pad=0.018,rounding_size=0.025",
        linewidth=1.4,
        edgecolor=edge,
        facecolor=face,
    )
    ax.add_patch(patch)
    ax.text(
        xy[0] + width / 2,
        xy[1] + height / 2,
        text,
        ha="center",
        va="center",
        color=NAVY,
        fontsize=fontsize,
        linespacing=1.25,
    )


def arrow(ax, start, end, color=GREY):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=12,
            linewidth=1.25,
            color=color,
        )
    )


def figure1() -> None:
    support = json.loads(
        (ROOT / "01_support_repair_benchmark_results.json").read_text(encoding="utf-8")
    )
    cases = support["cases"]
    exact = [
        {"case": r["case"], "added_cost": r["exact"]["added_cost"]}
        for r in cases
    ]
    greedy = [
        {"case": r["case"], "added_cost": r["greedy"]["added_cost"]}
        for r in cases
    ]

    fig = plt.figure(figsize=(13.2, 6.8), constrained_layout=True)
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1])
    ax = fig.add_subplot(gs[0, 0])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")
    ax.set_title("(a) Claim-adaptive decision loop", loc="left", fontweight="bold")

    box(
        ax,
        (0.08, 0.82),
        0.84,
        0.11,
        "Declare target contrast, minimum effect,\ncandidate costs, and uncertainty set",
        NAVY,
        LIGHT,
    )
    arrow(ax, (0.50, 0.82), (0.50, 0.74))
    box(
        ax,
        (0.22, 0.62),
        0.56,
        0.11,
        "Gate 1: contrast in design row space?",
        BLUE,
        "#EEF5FA",
    )
    arrow(ax, (0.22, 0.675), (0.07, 0.675))
    ax.text(0.145, 0.70, "No", ha="center", color=RED, fontweight="bold")
    box(
        ax,
        (0.01, 0.48),
        0.30,
        0.11,
        "Minimum-cost\nsupport repair",
        RED,
        "#FFF3EF",
    )
    arrow(ax, (0.16, 0.48), (0.16, 0.38))
    box(
        ax,
        (0.01, 0.24),
        0.30,
        0.11,
        "Acquire rows that span\nthe declared contrast",
        ORANGE,
        "#FFF8EE",
        8.9,
    )
    arrow(ax, (0.31, 0.295), (0.50, 0.62))

    arrow(ax, (0.78, 0.675), (0.92, 0.675))
    ax.text(0.85, 0.70, "Yes", ha="center", color=GREEN, fontweight="bold")
    box(
        ax,
        (0.69, 0.48),
        0.30,
        0.11,
        "Gate 2: interval excludes\nthe decision boundary?",
        GREEN,
        "#EEF8F4",
        8.9,
    )
    arrow(ax, (0.84, 0.48), (0.84, 0.38))
    box(
        ax,
        (0.69, 0.24),
        0.30,
        0.11,
        "Reduce contrast uncertainty:\nrefine, repeat, or audit",
        PURPLE,
        "#F5F1FA",
        8.8,
    )
    arrow(ax, (0.69, 0.295), (0.50, 0.62))

    box(
        ax,
        (0.28, 0.04),
        0.44,
        0.10,
        "Stop: effect present / effect absent /\nindeterminate",
        NAVY,
        "#EDF1F7",
        9.0,
    )
    arrow(ax, (0.84, 0.48), (0.61, 0.14), GREEN)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.set_title(
        "(b) Minimum-cost support repair across design families",
        loc="left",
        fontweight="bold",
    )
    labels = [
        r["case"]
        .replace("binary_", "")
        .replace("_", " ")
        .replace("alternative support", "alt. support")
        .replace("condition constrained", "condition")
        .replace("fractional pair contrasts", "multi-contrast")
        .replace("sparse two way model", "sparse 2-way")
        .replace("extreme level interaction", "extreme interaction")
        for r in exact
    ]
    exact_cost = np.array([float(r["added_cost"]) for r in exact])
    greedy_cost = np.array([float(r["added_cost"]) for r in greedy])
    y = np.arange(len(labels))
    h = 0.34
    ax2.barh(y + h / 2, exact_cost, height=h, color=BLUE, label="Exact")
    ax2.barh(y - h / 2, greedy_cost, height=h, color=ORANGE, label="Greedy")
    ax2.set_yticks(y, labels)
    ax2.invert_yaxis()
    ax2.set_xlabel("Added acquisition cost")
    ax2.grid(axis="x", color="#D6DCE3", linewidth=0.7)
    ax2.set_axisbelow(True)
    ax2.legend(frameon=False, ncol=2, loc="lower right")
    for i, (e, g) in enumerate(zip(exact_cost, greedy_cost)):
        if g > e + 1e-12:
            ax2.text(
                g + 0.35,
                i - h / 2,
                f"+{100 * (g / e - 1):.1f}%",
                va="center",
                color=RED,
                fontweight="bold",
            )
    ax2.text(
        0.02,
        0.02,
        "Feasible contrasts: exact 6/6; greedy 6/6\nExact-cost match: 5/6",
        transform=ax2.transAxes,
        ha="left",
        va="bottom",
        color=NAVY,
        bbox=dict(boxstyle="round,pad=0.35", fc="white", ec="#BCC7D2"),
    )
    fig.suptitle(
        "From claim declaration to minimum-cost evidence acquisition",
        fontsize=13,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "261_fig_claim_adaptive_method")


def figure2() -> None:
    gate = json.loads(
        (ROOT / "02_gate2_calibration_summary.json").read_text(encoding="utf-8")
    )
    loop = json.loads(
        (ROOT / "04b_multilevel_closed_loop_summary.json").read_text(
            encoding="utf-8"
        )
    )
    normal = gate["normal_covariance_primary"]
    rhos = sorted({float(r["correlation"]) for r in normal})
    false_present = [
        next(
            r["false_resolved_fraction"]
            for r in normal
            if float(r["correlation"]) == rho and float(r["true_delta"]) == 3.0
        )
        * 100
        for rho in rhos
    ]
    false_absent = [
        next(
            r["false_unresolved_fraction"]
            for r in normal
            if float(r["correlation"]) == rho and float(r["true_delta"]) == 5.0
        )
        * 100
        for rho in rhos
    ]

    rows = loop["rows"]
    policy_order = [
        "proposed_decision_targeted_switch",
        "bound_reduction_ablation",
        "c_optimal",
        "D_optimal",
        "prediction_variance",
        "high_fidelity_heavy",
        "random",
    ]
    short = {
        "proposed_decision_targeted_switch": "Proposed",
        "bound_reduction_ablation": "Bound ablation",
        "c_optimal": "c-optimal",
        "D_optimal": "D-optimal",
        "prediction_variance": "Prediction variance",
        "high_fidelity_heavy": "HF-heavy",
        "random": "Random",
    }

    fig, axes = plt.subplots(1, 3, figsize=(14.2, 4.9), constrained_layout=True)
    ax = axes[0]
    x = np.arange(len(rhos))
    w = 0.36
    ax.bar(x - w / 2, false_present, w, color=RED, label="False present, Δ=3")
    ax.bar(x + w / 2, false_absent, w, color=PURPLE, label="False absent, Δ=5")
    ax.axhline(5, color=GREY, linestyle="--", linewidth=1, label="5% reference")
    ax.set_xticks(x, [str(r) for r in rhos])
    ax.set_xlabel("Cell-error correlation ρ")
    ax.set_ylabel("Decision frequency (%)")
    ax.set_title("(a) Covariance-rule calibration", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8, loc="upper right")
    ax.grid(axis="y", color="#D6DCE3", linewidth=0.7)

    for j, delta in enumerate([2.0, 5.0], start=1):
        ax = axes[j]
        subset = {
            r["policy"]: r for r in rows if float(r["true_delta"]) == delta
        }
        costs = [subset[p]["median_cost_to_valid_decision"] for p in policy_order]
        valid = [subset[p]["valid_decision_fraction"] * 100 for p in policy_order]
        colors = [GREEN] + [CYAN, BLUE, BLUE, BLUE, ORANGE, GREY]
        y = np.arange(len(policy_order))
        ax.barh(y, costs, color=colors)
        ax.set_yticks(y, [short[p] for p in policy_order])
        ax.invert_yaxis()
        ax.set_xlabel("Median cost to valid decision")
        ax.set_title(
            f"({chr(97+j)}) True contrast Δ={int(delta)}",
            loc="left",
            fontweight="bold",
        )
        ax.grid(axis="x", color="#D6DCE3", linewidth=0.7)
        xmax = max(costs) * 1.32
        ax.set_xlim(0, xmax)
        for i, (cost, vf) in enumerate(zip(costs, valid)):
            ax.text(
                cost + 0.5,
                i,
                f"{vf:.1f}% valid",
                va="center",
                fontsize=8.3,
                color=NAVY,
            )
        if delta == 2:
            ax.annotate(
                "35% lower than\nvariance baselines",
                xy=(19.5, 0),
                xytext=(25, 1.2),
                arrowprops=dict(arrowstyle="->", color=GREEN),
                color=GREEN,
                fontweight="bold",
                fontsize=8.5,
            )

    fig.suptitle(
        "Calibrated stopping and cost to a valid decision",
        fontsize=13,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "262_fig_cost_to_valid_decision")


def figure3() -> None:
    data = json.loads(
        (ROOT / "07_joint_mapping_ensemble_audit.json").read_text(encoding="utf-8")
    )["results"]
    labels = [
        f"{r['case_id'].replace('alph15-', '')}\n"
        f"{'fine' if r['support'] == 'fine_F100' else 'shared'}"
        for r in data
    ]
    joint = [r["ensemble_resolution_ratio"] for r in data]
    boxv = [r["axis_aligned_resolution_ratio"] for r in data]
    contraction = [100 * r["joint_to_axis_aligned_diameter_ratio"] for r in data]
    x = np.arange(len(labels))

    fig, axes = plt.subplots(
        1, 2, figsize=(11.8, 4.8), gridspec_kw={"width_ratios": [1.35, 0.8]},
        constrained_layout=True
    )
    ax = axes[0]
    w = 0.36
    ax.bar(x - w / 2, joint, w, color=GREEN, label="Joint pipeline ensemble")
    ax.bar(x + w / 2, boxv, w, color=PURPLE, label="Independent-cell box")
    ax.axhline(3, color=RED, linestyle="--", linewidth=1.4, label="Frozen threshold")
    ax.set_xticks(x, labels)
    ax.set_ylabel("Mapping-robustness ratio")
    ax.set_title("(a) Opposite decisions from the same mapping outputs", loc="left", fontweight="bold")
    ax.legend(frameon=False, fontsize=8.5)
    ax.grid(axis="y", color="#D6DCE3", linewidth=0.7)
    for i, value in enumerate(joint):
        ax.text(i - w / 2, value + 0.16, f"{value:.2f}", ha="center", fontsize=8.4)
    for i, value in enumerate(boxv):
        ax.text(i + w / 2, value + 0.16, f"{value:.2f}", ha="center", fontsize=8.4)

    ax2 = axes[1]
    ax2.bar(x, contraction, color=ORANGE)
    ax2.set_xticks(x, labels)
    ax2.set_ylabel("Joint diameter / box envelope (%)")
    ax2.set_ylim(0, 10)
    ax2.set_title("(b) Dependence-preserving contraction", loc="left", fontweight="bold")
    ax2.grid(axis="y", color="#D6DCE3", linewidth=0.7)
    for i, value in enumerate(contraction):
        ax2.text(i, value + 0.2, f"{value:.1f}%", ha="center", fontweight="bold", color=NAVY)

    fig.suptitle(
        "Uncertainty-set geometry reversed the real-CFD mapping decision",
        fontsize=13,
        fontweight="bold",
        color=NAVY,
    )
    save(fig, "263_fig_joint_vs_box_mapping")


if __name__ == "__main__":
    figure1()
    figure2()
    figure3()
    print("Built V9 main figures 261-263 (PNG and PDF).")
