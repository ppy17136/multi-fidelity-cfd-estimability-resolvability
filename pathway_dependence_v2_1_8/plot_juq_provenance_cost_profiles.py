"""Generate the provenance cost-profile figure from the released data."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np


HERE = Path(__file__).resolve().parent
INPUT = HERE / "37_REAL_CFD_COMMON_PARTITION_PROFILES.json"
PNG = HERE / "figures" / "43_FIG_PROVENANCE_COST_PROFILES.png"
PDF = HERE / "figures" / "43_FIG_PROVENANCE_COST_PROFILES.pdf"
SVG = HERE / "figures" / "43_FIG_PROVENANCE_COST_PROFILES.svg"

GEOMETRIES = {
    "alph05-10071-2024": "G1",
    "alph15-10929-3036": "G2",
    "alph15-13929-2024": "G3",
    "alph15-7929-2024": "G4",
}
SUPPORT = {"fine_F100": "fine", "shared_F000": "shared"}
EPS_KEYS = ("0e+00", "1e-10", "1e-09", "1e-08", "1e-07")
EPS_LABELS = ("exact", r"$10^{-10}$", r"$10^{-9}$", r"$10^{-8}$", r"$10^{-7}$")


def model(audit, name):
    return next(item for item in audit["cost_models"] if item["cost_model"] == name)


def main():
    payload = json.loads(INPUT.read_text(encoding="utf-8"))
    audits = payload["results"]

    plt.rcParams.update({
        "font.family": "DejaVu Sans",
        "font.size": 10.0,
        "axes.titlesize": 10.5,
        "axes.labelsize": 9.5,
        "legend.fontsize": 9.0,
        "xtick.labelsize": 8.5,
        "ytick.labelsize": 8.5,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
    })

    fig = plt.figure(figsize=(12.0, 7.6), constrained_layout=True)
    grid = fig.add_gridspec(2, 2, width_ratios=(0.86, 1.7), height_ratios=(1, 1))
    ax_graph = fig.add_subplot(grid[0, 0])
    ax_table = fig.add_subplot(grid[1, 0])
    ax_ratio = fig.add_subplot(grid[0, 1])
    ax_heat = fig.add_subplot(grid[1, 1])

    # (a) Declared mesh-by-closure provenance graph.
    positions = {
        "F000": (0.15, 0.78), "F100": (0.85, 0.78),
        "F010": (0.15, 0.22), "F110": (0.85, 0.22),
    }
    mesh_color, closure_color = "#2667A8", "#D97706"
    for left, right in (("F000", "F100"), ("F010", "F110")):
        x = [positions[left][0], positions[right][0]]
        y = [positions[left][1], positions[right][1]]
        ax_graph.plot(x, y, color=mesh_color, lw=3.0, solid_capstyle="round")
    for left, right in (("F000", "F010"), ("F100", "F110")):
        x = [positions[left][0], positions[right][0]]
        y = [positions[left][1], positions[right][1]]
        ax_graph.plot(x, y, color=closure_color, lw=3.0, solid_capstyle="round")
    for name, (x, y) in positions.items():
        ax_graph.scatter(x, y, s=780, color="white", edgecolor="#222222", linewidth=1.4, zorder=3)
        ax_graph.text(x, y, name, ha="center", va="center", weight="bold", zorder=4)
    ax_graph.text(0.50, 0.86, "mesh", color=mesh_color, ha="center", weight="bold")
    ax_graph.text(0.50, 0.14, "mesh", color=mesh_color, ha="center", weight="bold")
    ax_graph.text(0.06, 0.50, "closure", color=closure_color, ha="center", va="center", rotation=90, weight="bold")
    ax_graph.text(0.94, 0.50, "closure", color=closure_color, ha="center", va="center", rotation=90, weight="bold")
    ax_graph.text(0.50, 0.02, "Cut edges permit pathway disagreement", ha="center", color="#444444")
    ax_graph.set_xlim(-0.05, 1.05)
    ax_graph.set_ylim(-0.04, 1.03)
    ax_graph.axis("off")
    ax_graph.set_title("(a) Declared provenance graph", loc="left", weight="bold")

    # Cost-model table kept separate to prevent labels crowding the graph.
    ax_table.axis("off")
    cell_text = [
        ["Balanced", "0.25", "0.25"],
        ["Mesh-dominant", "0.40", "0.10"],
        ["Closure-dominant", "0.10", "0.40"],
    ]
    table = ax_table.table(
        cellText=cell_text,
        colLabels=["Normalized edge cost", "mesh", "closure"],
        loc="upper center",
        cellLoc="center",
        colWidths=[0.54, 0.23, 0.23],
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8.6)
    table.scale(1.0, 1.45)
    for (row, col), cell in table.get_celld().items():
        cell.set_edgecolor("#B9C1CA")
        if row == 0:
            cell.set_facecolor("#E8EDF3")
            cell.set_text_props(weight="bold")
        elif row % 2 == 0:
            cell.set_facecolor("#F6F7F9")
    ax_table.text(
        0.0, 0.30,
        "Weights are declared sensitivities, not fitted probabilities\nor physical error magnitudes.",
        transform=ax_table.transAxes, ha="left", va="top", color="#444444",
    )
    ax_table.set_title("Cost models (total graph cost = 1)", loc="left", weight="bold")

    # (b) Common-partition ratio profiles under balanced weights.
    colors = plt.get_cmap("tab10").colors[:4]
    for audit in audits:
        geometry = GEOMETRIES[audit["case_id"]]
        geometry_index = int(geometry[1]) - 1
        support = SUPPORT[audit["support"]]
        profile = model(audit, "balanced")["profile"]
        x = np.asarray([row["relaxation_budget"] for row in profile])
        y = np.asarray([row["minimum_single_partition_ratio"] for row in profile])
        ax_ratio.step(
            x, y, where="post", color=colors[geometry_index],
            ls="-" if support == "fine" else "--",
            lw=1.9, marker="o", ms=3.7,
            label=f"{geometry}-{support}",
        )
    lower, upper = 1.6442596248017918, 2.2831676147449826
    ax_ratio.axhspan(lower, upper, color="#A7D8A6", alpha=0.30, zorder=0)
    ax_ratio.axvline(0.5, color="#4B5563", lw=1.0, ls=":")
    ax_ratio.text(0.505, 16.0, "$c=0.5$", rotation=90, va="top", ha="left", color="#4B5563")
    ax_ratio.text(0.98, np.sqrt(lower * upper), "common threshold interval", ha="right", va="center", color="#2F6B36")
    ax_ratio.set_yscale("log")
    ax_ratio.set_xlim(-0.02, 1.02)
    ax_ratio.set_ylim(0.11, 30)
    ax_ratio.set_xlabel("Normalized relaxation budget, $c$")
    ax_ratio.set_ylabel(r"Worst ratio, $R_*(c)$", labelpad=4)
    ax_ratio.grid(True, which="both", color="#D9DEE5", lw=0.6, alpha=0.75)
    ax_ratio.legend(ncol=4, loc="lower left", frameon=True, columnspacing=1.0, handlelength=2.2)
    ax_ratio.set_title("(b) Cost-indexed decision fragility across all eight audits", loc="left", weight="bold")

    # (c) Exact-to-tolerance diameter cost, balanced model.
    heat = np.asarray([
        [float(model(audit, "balanced")["minimum_diameter_cost_by_epsilon"][key]) for key in EPS_KEYS]
        for audit in audits
    ])
    image = ax_heat.imshow(heat, aspect="auto", cmap="viridis", norm=Normalize(0.0, 1.0))
    labels = [f"{GEOMETRIES[a['case_id']]}-{SUPPORT[a['support']]}" for a in audits]
    ax_heat.set_xticks(range(len(EPS_LABELS)), EPS_LABELS)
    ax_heat.set_yticks(range(len(labels)), labels)
    ax_heat.set_xlabel(r"Relative diameter tolerance, $\varepsilon$")
    ax_heat.set_ylabel("Audit", labelpad=4)
    for row in range(heat.shape[0]):
        for col in range(heat.shape[1]):
            value = heat[row, col]
            color = "white" if value < 0.45 or value > 0.78 else "black"
            ax_heat.text(col, row, f"{value:.1f}", ha="center", va="center", color=color, weight="bold", fontsize=8.0)
    bar = fig.colorbar(image, ax=ax_heat, fraction=0.04, pad=0.02)
    bar.set_label(r"Cost $\rho_D(\varepsilon)$", labelpad=5)
    ax_heat.set_title("(c) Exact equality is less stable than near-attainment", loc="left", weight="bold")

    fig.suptitle(
        "Cost-indexed provenance relaxation in multifidelity UQ",
        fontsize=12.5, weight="bold",
    )
    fig.savefig(PNG, dpi=360, bbox_inches="tight", facecolor="white")
    fig.savefig(PDF, bbox_inches="tight", facecolor="white")
    fig.savefig(SVG, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    print(PNG)
    print(PDF)
    print(SVG)


if __name__ == "__main__":
    main()
