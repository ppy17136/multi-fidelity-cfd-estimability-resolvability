#!/usr/bin/env python3
"""Manufactured validation of the two-gate estimability-resolvability rule."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import numpy as np
from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
SEED = 20260728
N_MONTE_CARLO = 200_000
CELL_ERROR_BOUND = 0.25
CONTRAST_FLOOR = 4.0 * CELL_ERROR_BOUND
MATERIALITY_THRESHOLD = 3.0
RANK_TOLERANCE = 1.0e-10

CSV_PATH = ROOT / "199_manufactured_two_gate_summary.csv"
JSON_PATH = ROOT / "199_manufactured_two_gate_summary.json"
QUANTILE_CSV_PATH = ROOT / "199_manufactured_two_gate_quantiles.csv"
PNG_PATH = ROOT / "199_fig_manufactured_two_gate_benchmark.png"
SVG_PATH = ROOT / "199_fig_manufactured_two_gate_benchmark.svg"
OK_PATH = ROOT / "MANUFACTURED_TWO_GATE_BENCHMARK_199_COMPLETE.ok"

CELL_ORDER = ["00", "10", "01", "11"]
FULL_X = np.array(
    [
        [1.0, 0.0, 0.0, 0.0],
        [1.0, 1.0, 0.0, 0.0],
        [1.0, 0.0, 1.0, 0.0],
        [1.0, 1.0, 1.0, 1.0],
    ]
)
MISSING_11_X = FULL_X[:3, :]
COEFFICIENT_CONTRAST = np.array([0.0, 0.0, 0.0, 1.0])
CELL_CONTRAST = np.array([1.0, -1.0, -1.0, 1.0])


def estimability_audit(x: np.ndarray) -> dict[str, object]:
    rank = int(np.linalg.matrix_rank(x))
    projected = np.linalg.pinv(x) @ x @ COEFFICIENT_CONTRAST
    residual = COEFFICIENT_CONTRAST - projected
    relative_residual = float(
        np.linalg.norm(residual)
        / max(float(np.linalg.norm(COEFFICIENT_CONTRAST)), 1.0)
    )
    estimable = relative_residual <= RANK_TOLERANCE

    weights: list[float] | None = None
    weight_residual: float | None = None
    if estimable:
        solved, *_ = np.linalg.lstsq(x.T, COEFFICIENT_CONTRAST, rcond=None)
        weights = [float(value) for value in solved]
        weight_residual = float(
            np.linalg.norm(solved @ x - COEFFICIENT_CONTRAST)
        )

    return {
        "rank": rank,
        "row_space_relative_residual": relative_residual,
        "estimable": estimable,
        "estimator_weights": weights,
        "estimator_weight_residual": weight_residual,
    }


def ratio_quantiles(values: np.ndarray) -> dict[str, float]:
    probabilities = [0.0, 0.01, 0.05, 0.25, 0.50, 0.75, 0.95, 0.99, 1.0]
    computed = np.quantile(values, probabilities)
    return {
        f"q{int(probability * 100):02d}": float(value)
        for probability, value in zip(probabilities, computed, strict=True)
    }


def run_complete_scenario(
    rng: np.random.Generator,
    name: str,
    true_interaction: float,
    expected_decision: str,
) -> tuple[dict[str, object], np.ndarray]:
    beta = np.array([10.0, 1.5, -0.5, true_interaction])
    true_cell_response = FULL_X @ beta

    perturbations = rng.uniform(
        -CELL_ERROR_BOUND,
        CELL_ERROR_BOUND,
        size=(N_MONTE_CARLO, 4),
    )
    contrast_errors = perturbations @ CELL_CONTRAST
    observed_contrasts = true_interaction + contrast_errors
    ratios = np.abs(observed_contrasts) / CONTRAST_FLOOR
    passes = ratios > MATERIALITY_THRESHOLD

    audit = estimability_audit(FULL_X)
    result: dict[str, object] = {
        "scenario": name,
        "support": "complete_2x2",
        **audit,
        "true_interaction": true_interaction,
        "true_ratio": abs(true_interaction) / CONTRAST_FLOOR,
        "expected_decision": expected_decision,
        "numerical_floor": CONTRAST_FLOOR,
        "cell_error_bound": CELL_ERROR_BOUND,
        "threshold": MATERIALITY_THRESHOLD,
        "monte_carlo_repeats": N_MONTE_CARLO,
        "maximum_observed_abs_contrast_error": float(
            np.max(np.abs(contrast_errors))
        ),
        "mean_observed_ratio": float(np.mean(ratios)),
        "pass_rate": float(np.mean(passes)),
        "nonpass_rate": float(np.mean(~passes)),
        "ratio_quantiles": ratio_quantiles(ratios),
        "true_cell_response": {
            cell: float(value)
            for cell, value in zip(CELL_ORDER, true_cell_response, strict=True)
        },
    }
    return result, ratios


def write_csv(rows: list[dict[str, object]]) -> None:
    fieldnames = [
        "scenario",
        "support",
        "rank",
        "row_space_relative_residual",
        "estimable",
        "true_interaction",
        "true_ratio",
        "expected_decision",
        "numerical_floor",
        "cell_error_bound",
        "threshold",
        "monte_carlo_repeats",
        "maximum_observed_abs_contrast_error",
        "mean_observed_ratio",
        "pass_rate",
        "nonpass_rate",
    ]
    with CSV_PATH.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(stream, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({key: row.get(key) for key in fieldnames})


def write_quantile_csv(rows: list[dict[str, object]]) -> None:
    quantile_names = ["q00", "q01", "q05", "q25", "q50", "q75", "q95", "q99", "q100"]
    with QUANTILE_CSV_PATH.open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.DictWriter(
            stream,
            fieldnames=["scenario", "true_ratio", *quantile_names],
        )
        writer.writeheader()
        for row in rows:
            if not row.get("ratio_quantiles"):
                continue
            writer.writerow(
                {
                    "scenario": row["scenario"],
                    "true_ratio": row["true_ratio"],
                    **row["ratio_quantiles"],
                }
            )


def load_font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    family = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    try:
        return ImageFont.truetype(family, size=size)
    except OSError:
        pass
    return ImageFont.load_default()


def draw_lattice(
    draw: ImageDraw.ImageDraw,
    origin: tuple[int, int],
    missing_11: bool,
    title: str,
    decision: str,
) -> None:
    x0, y0 = origin
    size = 155
    draw.rectangle((x0, y0, x0 + size, y0 + size), outline="#4c566a", width=4)
    coordinates = {
        "00": (x0, y0 + size),
        "10": (x0 + size, y0 + size),
        "01": (x0, y0),
        "11": (x0 + size, y0),
    }
    for cell, (x, y) in coordinates.items():
        missing = missing_11 and cell == "11"
        fill = "#ffffff" if missing else "#4c78a8"
        outline = "#c44e52" if missing else "#1f3b5c"
        draw.ellipse((x - 13, y - 13, x + 13, y + 13), fill=fill, outline=outline, width=4)
        draw.text((x + 17, y - 12), cell, font=load_font(21), fill="#202531")
        if missing:
            draw.line((x - 9, y - 9, x + 9, y + 9), fill="#c44e52", width=4)
            draw.line((x - 9, y + 9, x + 9, y - 9), fill="#c44e52", width=4)
    draw.text((x0 - 5, y0 - 62), title, font=load_font(25, bold=True), fill="#202531")
    color = "#b23a48" if "FAIL" in decision else "#237a57"
    draw.text((x0 - 5, y0 + size + 35), decision, font=load_font(24, bold=True), fill=color)


def build_png(
    low: dict[str, object],
    boundary: dict[str, object],
    high: dict[str, object],
) -> None:
    width, height = 1800, 1050
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    draw.text(
        (60, 35),
        "Manufactured validation of the two-gate estimability-resolvability rule",
        font=load_font(40, bold=True),
        fill="#1f2430",
    )

    draw.text((60, 112), "A  Gate 1: support and estimability", font=load_font(31, bold=True), fill="#1f2430")
    draw_lattice(draw, (115, 245), True, "Missing 11 corner", "Gate 1 FAIL")
    draw_lattice(draw, (430, 245), False, "Complete 2 x 2", "Gate 1 PASS")
    draw.text(
        (70, 520),
        "Missing support cannot identify the AB interaction;\n"
        "an emulator prediction at 11 does not repair estimability.",
        font=load_font(24),
        fill="#424854",
        spacing=8,
    )

    draw.text((760, 112), "B  Gate 2: numerical resolvability", font=load_font(31, bold=True), fill="#1f2430")
    chart_left, chart_top = 790, 225
    chart_right, chart_bottom = 1700, 700
    x_min, x_max = 0.0, 6.5

    def x_pos(value: float) -> int:
        return int(chart_left + (value - x_min) / (x_max - x_min) * (chart_right - chart_left))

    draw.line((chart_left, chart_bottom, chart_right, chart_bottom), fill="#333333", width=3)
    for tick in range(0, 7):
        x = x_pos(float(tick))
        draw.line((x, chart_bottom, x, chart_bottom + 10), fill="#333333", width=2)
        draw.text((x - 7, chart_bottom + 18), str(tick), font=load_font(20), fill="#333333")
    draw.text(
        (1110, chart_bottom + 65),
        "Observed contrast / numerical evidence floor",
        font=load_font(25, bold=True),
        fill="#333333",
    )

    threshold_x = x_pos(MATERIALITY_THRESHOLD)
    draw.line((threshold_x, chart_top - 20, threshold_x, chart_bottom), fill="#c44e52", width=5)
    draw.text(
        (threshold_x + 10, chart_top - 22),
        "threshold = 3",
        font=load_font(22, bold=True),
        fill="#c44e52",
    )

    scenario_rows = [
        (low, 320, "#4c78a8", "Below-floor control"),
        (boundary, 455, "#f2a541", "Threshold sensitivity"),
        (high, 590, "#54a66f", "Above-floor control"),
    ]
    for row, y, color, label in scenario_rows:
        q = row["ratio_quantiles"]
        q00, q05, q50, q95, q100 = (
            float(q["q00"]),
            float(q["q05"]),
            float(q["q50"]),
            float(q["q95"]),
            float(q["q100"]),
        )
        draw.text((chart_left, y - 62), label, font=load_font(23, bold=True), fill="#282d38")
        draw.line((x_pos(q00), y, x_pos(q100), y), fill=color, width=8)
        draw.line((x_pos(q05), y, x_pos(q95), y), fill=color, width=22)
        draw.ellipse(
            (x_pos(q50) - 9, y - 9, x_pos(q50) + 9, y + 9),
            fill="#ffffff",
            outline="#202531",
            width=4,
        )
        draw.text(
            (x_pos(q100) + 14, y - 18),
            f"pass rate {float(row['pass_rate']) * 100:.1f}%",
            font=load_font(21),
            fill="#282d38",
        )

    draw.text((60, 825), "C  Predeclared decision logic", font=load_font(31, bold=True), fill="#1f2430")
    decision_text = (
        "Missing corner: do not estimate or interpret the interaction.    "
        "Complete + 2F: estimable but not resolved.    "
        "Complete + 5F: estimable and resolved."
    )
    draw.rounded_rectangle((60, 885, 1740, 995), radius=22, fill="#f1f3f7", outline="#8792a2", width=3)
    draw.text((90, 920), decision_text, font=load_font(25, bold=True), fill="#253147")

    image.save(PNG_PATH, dpi=(300, 300))


def build_svg(
    low: dict[str, object],
    boundary: dict[str, object],
    high: dict[str, object],
) -> None:
    def esc(text: str) -> str:
        return (
            text.replace("&", "&amp;")
            .replace("<", "&lt;")
            .replace(">", "&gt;")
            .replace('"', "&quot;")
        )

    x_left, x_right = 760.0, 1660.0

    def x_pos(value: float) -> float:
        return x_left + value / 6.5 * (x_right - x_left)

    parts = [
        '<svg xmlns="http://www.w3.org/2000/svg" width="1800" height="1050" viewBox="0 0 1800 1050">',
        '<rect width="1800" height="1050" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#1f2430}.title{font-size:40px;font-weight:700}.head{font-size:31px;font-weight:700}.label{font-size:24px}.bold{font-weight:700}</style>',
        '<text x="60" y="70" class="title">Manufactured validation of the two-gate estimability-resolvability rule</text>',
        '<text x="60" y="145" class="head">A  Gate 1: support and estimability</text>',
    ]

    for origin_x, missing, title, decision in [
        (115, True, "Missing 11 corner", "Gate 1 FAIL"),
        (430, False, "Complete 2 x 2", "Gate 1 PASS"),
    ]:
        y0, size = 245, 155
        parts.append(f'<text x="{origin_x}" y="205" class="label bold">{esc(title)}</text>')
        parts.append(f'<rect x="{origin_x}" y="{y0}" width="{size}" height="{size}" fill="none" stroke="#4c566a" stroke-width="4"/>')
        points = {
            "00": (origin_x, y0 + size),
            "10": (origin_x + size, y0 + size),
            "01": (origin_x, y0),
            "11": (origin_x + size, y0),
        }
        for cell, (x, y) in points.items():
            is_missing = missing and cell == "11"
            fill = "white" if is_missing else "#4c78a8"
            stroke = "#c44e52" if is_missing else "#1f3b5c"
            parts.append(f'<circle cx="{x}" cy="{y}" r="13" fill="{fill}" stroke="{stroke}" stroke-width="4"/>')
            parts.append(f'<text x="{x + 18}" y="{y + 7}" class="label">{cell}</text>')
            if is_missing:
                parts.append(f'<path d="M{x-9},{y-9} L{x+9},{y+9} M{x-9},{y+9} L{x+9},{y-9}" stroke="#c44e52" stroke-width="4"/>')
        decision_color = "#b23a48" if "FAIL" in decision else "#237a57"
        parts.append(f'<text x="{origin_x}" y="465" class="label bold" style="fill:{decision_color}">{decision}</text>')

    parts.extend(
        [
            '<text x="70" y="545" class="label">Missing support cannot identify the AB interaction;</text>',
            '<text x="70" y="580" class="label">an emulator prediction at 11 does not repair estimability.</text>',
            '<text x="730" y="145" class="head">B  Gate 2: numerical resolvability</text>',
            f'<line x1="{x_left}" y1="700" x2="{x_right}" y2="700" stroke="#333" stroke-width="3"/>',
        ]
    )
    for tick in range(7):
        x = x_pos(float(tick))
        parts.append(f'<line x1="{x}" y1="700" x2="{x}" y2="710" stroke="#333" stroke-width="2"/>')
        parts.append(f'<text x="{x - 6}" y="738" class="label">{tick}</text>')
    parts.append('<text x="1050" y="790" class="label bold">Observed contrast / numerical evidence floor</text>')
    tx = x_pos(MATERIALITY_THRESHOLD)
    parts.append(f'<line x1="{tx}" y1="205" x2="{tx}" y2="700" stroke="#c44e52" stroke-width="5"/>')
    parts.append(f'<text x="{tx + 10}" y="225" class="label bold" style="fill:#c44e52">threshold = 3</text>')

    for row, y, color, label in [
        (low, 320, "#4c78a8", "Below-floor control"),
        (boundary, 455, "#f2a541", "Threshold sensitivity"),
        (high, 590, "#54a66f", "Above-floor control"),
    ]:
        q = row["ratio_quantiles"]
        parts.append(f'<text x="{x_left}" y="{y - 40}" class="label bold">{esc(label)}</text>')
        parts.append(f'<line x1="{x_pos(float(q["q00"]))}" y1="{y}" x2="{x_pos(float(q["q100"]))}" y2="{y}" stroke="{color}" stroke-width="8"/>')
        parts.append(f'<line x1="{x_pos(float(q["q05"]))}" y1="{y}" x2="{x_pos(float(q["q95"]))}" y2="{y}" stroke="{color}" stroke-width="22"/>')
        parts.append(f'<circle cx="{x_pos(float(q["q50"]))}" cy="{y}" r="9" fill="white" stroke="#202531" stroke-width="4"/>')
        parts.append(f'<text x="{x_pos(float(q["q100"])) + 14}" y="{y + 7}" class="label">pass rate {float(row["pass_rate"]) * 100:.1f}%</text>')

    parts.extend(
        [
            '<text x="60" y="855" class="head">C  Predeclared decision logic</text>',
            '<rect x="60" y="885" width="1680" height="110" rx="22" fill="#f1f3f7" stroke="#8792a2" stroke-width="3"/>',
            '<text x="90" y="950" class="label bold">Missing corner: do not estimate.   Complete + 2F: estimable, not resolved.   Complete + 5F: estimable and resolved.</text>',
            '</svg>',
        ]
    )
    SVG_PATH.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    rng = np.random.default_rng(SEED)

    missing_audit = estimability_audit(MISSING_11_X)
    missing_result: dict[str, object] = {
        "scenario": "missing_corner_high_signal",
        "support": "missing_11",
        **missing_audit,
        "true_interaction": 5.0,
        "true_ratio": 5.0,
        "expected_decision": "gate1_fail_gate2_not_evaluated",
        "numerical_floor": CONTRAST_FLOOR,
        "cell_error_bound": CELL_ERROR_BOUND,
        "threshold": MATERIALITY_THRESHOLD,
        "monte_carlo_repeats": 0,
        "maximum_observed_abs_contrast_error": None,
        "mean_observed_ratio": None,
        "pass_rate": None,
        "nonpass_rate": None,
        "ratio_quantiles": None,
    }

    low_result, low_ratios = run_complete_scenario(
        rng,
        "complete_below_floor",
        true_interaction=2.0,
        expected_decision="gate1_pass_gate2_fail",
    )
    boundary_result, boundary_ratios = run_complete_scenario(
        rng,
        "complete_at_threshold",
        true_interaction=3.0,
        expected_decision="transition_region",
    )
    high_result, high_ratios = run_complete_scenario(
        rng,
        "complete_above_floor",
        true_interaction=5.0,
        expected_decision="gate1_pass_gate2_pass",
    )

    rows = [missing_result, low_result, boundary_result, high_result]

    assert missing_result["rank"] == 3
    assert not missing_result["estimable"]
    assert low_result["rank"] == 4 and low_result["estimable"]
    assert high_result["rank"] == 4 and high_result["estimable"]
    assert float(low_result["pass_rate"]) == 0.0
    assert float(high_result["pass_rate"]) == 1.0
    assert abs(float(boundary_result["pass_rate"]) - 0.5) < 0.01
    assert np.max(np.abs(low_ratios - 2.0)) <= 1.0 + 1.0e-12
    assert np.max(np.abs(high_ratios - 5.0)) <= 1.0 + 1.0e-12
    assert np.max(np.abs(boundary_ratios - 3.0)) <= 1.0 + 1.0e-12

    write_csv(rows)
    write_quantile_csv(rows)
    JSON_PATH.write_text(
        json.dumps(
            {
                "seed": SEED,
                "monte_carlo_repeats": N_MONTE_CARLO,
                "cell_error_bound": CELL_ERROR_BOUND,
                "contrast_numerical_floor": CONTRAST_FLOOR,
                "materiality_threshold": MATERIALITY_THRESHOLD,
                "rank_tolerance": RANK_TOLERANCE,
                "cell_order": CELL_ORDER,
                "coefficient_contrast": COEFFICIENT_CONTRAST.tolist(),
                "cell_contrast": CELL_CONTRAST.tolist(),
                "scenarios": rows,
            },
            indent=2,
        ),
        encoding="utf-8",
    )
    build_png(low_result, boundary_result, high_result)
    build_svg(low_result, boundary_result, high_result)
    OK_PATH.write_text(
        "manufactured two-gate benchmark completed\n"
        f"seed={SEED}\n"
        f"monte_carlo_repeats={N_MONTE_CARLO}\n",
        encoding="utf-8",
    )

    print(f"wrote {CSV_PATH}")
    print(f"wrote {QUANTILE_CSV_PATH}")
    print(f"wrote {JSON_PATH}")
    print(f"wrote {PNG_PATH}")
    print(f"wrote {SVG_PATH}")
    for row in rows:
        print(
            row["scenario"],
            f"rank={row['rank']}",
            f"estimable={row['estimable']}",
            f"pass_rate={row['pass_rate']}",
        )


if __name__ == "__main__":
    main()
