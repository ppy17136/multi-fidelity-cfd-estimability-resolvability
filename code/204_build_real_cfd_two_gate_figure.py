#!/usr/bin/env python3
"""Build the real-CFD two-gate result figure."""

from pathlib import Path

from PIL import Image, ImageDraw, ImageFont


ROOT = Path(__file__).resolve().parent
PNG_PATH = ROOT / "204_fig_real_cfd_estimable_but_unresolved.png"
SVG_PATH = ROOT / "204_fig_real_cfd_estimable_but_unresolved.svg"

THRESHOLD = 3.0
ROWS = [
    {
        "geometry": "alph05-10071-2024",
        "contrast": 0.0057742,
        "floor": 0.0031035,
        "ratio": 1.861,
    },
    {
        "geometry": "alph15-10929-3036",
        "contrast": 0.0065753,
        "floor": 0.0033278,
        "ratio": 1.976,
    },
    {
        "geometry": "alph15-13929-2024",
        "contrast": 0.0062830,
        "floor": 0.0022779,
        "ratio": 2.758,
    },
    {
        "geometry": "alph15-7929-2024",
        "contrast": 0.0089263,
        "floor": 0.0033618,
        "ratio": 2.655,
    },
]
POOLED = [
    ("Pooled archived floor", 2.232394),
    ("Pooled conservative floor", 1.913890),
]


def font(size: int, bold: bool = False) -> ImageFont.FreeTypeFont | ImageFont.ImageFont:
    family = "DejaVuSans-Bold.ttf" if bold else "DejaVuSans.ttf"
    try:
        return ImageFont.truetype(family, size=size)
    except OSError:
        pass
    return ImageFont.load_default()


def build_png() -> None:
    width, height = 1700, 1050
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image)

    draw.text(
        (60, 40),
        "Four complete CFD lattices are estimable but numerically unresolved",
        font=font(40, bold=True),
        fill="#202531",
    )
    draw.text(
        (60, 98),
        "Gate 1: PASS for all four contrasts    |    Gate 2 threshold: ratio > 3",
        font=font(26, bold=True),
        fill="#3c4858",
    )

    left, right = 510, 1550
    top, bottom = 215, 760
    x_max = 3.4

    def xpos(value: float) -> int:
        return int(left + value / x_max * (right - left))

    for tick in [0, 1, 2, 3]:
        x = xpos(float(tick))
        draw.line((x, top - 15, x, bottom), fill="#d8dde5", width=2)
        draw.text((x - 8, bottom + 18), str(tick), font=font(23), fill="#3a404b")

    threshold_x = xpos(THRESHOLD)
    draw.line((threshold_x, top - 25, threshold_x, bottom), fill="#c44e52", width=6)
    draw.text(
        (threshold_x - 80, top - 70),
        "threshold = 3",
        font=font(23, bold=True),
        fill="#c44e52",
    )

    bar_height = 64
    row_gap = 118
    colors = ["#4c78a8", "#4c78a8", "#7a5dc7", "#7a5dc7"]
    for index, row in enumerate(ROWS):
        y = top + index * row_gap
        draw.text(
            (60, y + 15),
            row["geometry"],
            font=font(25, bold=True),
            fill="#202531",
        )
        draw.text(
            (60, y + 52),
            f"contrast {row['contrast']:.7f}  |  floor {row['floor']:.7f}",
            font=font(20),
            fill="#596273",
        )
        draw.rounded_rectangle(
            (left, y, xpos(row["ratio"]), y + bar_height),
            radius=14,
            fill=colors[index],
        )
        draw.text(
            (xpos(row["ratio"]) + 18, y + 13),
            f"{row['ratio']:.3f}",
            font=font(26, bold=True),
            fill="#202531",
        )
        draw.text(
            (right - 125, y + 13),
            "UNRESOLVED",
            font=font(21, bold=True),
            fill="#a13f48",
        )

    draw.line((60, 810, 1640, 810), fill="#aeb6c2", width=3)
    draw.text((60, 840), "Pooled checks", font=font(27, bold=True), fill="#202531")

    pooled_x0 = 420
    for index, (label, value) in enumerate(POOLED):
        y = 900 + index * 55
        draw.text((60, y - 14), label, font=font(21), fill="#3a404b")
        x = pooled_x0 + index * 530
        draw.polygon(
            [(x, y), (x + 13, y - 13), (x + 26, y), (x + 13, y + 13)],
            fill="#4c78a8" if index == 0 else "#7a5dc7",
        )
        draw.text((x + 42, y - 17), f"ratio {value:.3f} < 3", font=font(23, bold=True), fill="#3a404b")

    draw.text(
        (1030, 922),
        "Below threshold means insufficient audited evidence,\nnot a zero physical interaction.",
        font=font(21, bold=True),
        fill="#596273",
        spacing=7,
    )

    image.save(PNG_PATH, dpi=(300, 300))


def build_svg() -> None:
    width, height = 1700, 1050
    left, right = 510.0, 1550.0
    top, bottom = 215.0, 760.0
    x_max = 3.4

    def xpos(value: float) -> float:
        return left + value / x_max * (right - left)

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}" viewBox="0 0 {width} {height}">',
        '<rect width="1700" height="1050" fill="white"/>',
        '<style>text{font-family:Arial,sans-serif;fill:#202531}.title{font-size:40px;font-weight:700}.sub{font-size:26px;font-weight:700}.label{font-size:25px}.small{font-size:20px}.bold{font-weight:700}</style>',
        '<text x="60" y="75" class="title">Four complete CFD lattices are estimable but numerically unresolved</text>',
        '<text x="60" y="130" class="sub">Gate 1: PASS for all four contrasts | Gate 2 threshold: ratio &gt; 3</text>',
    ]

    for tick in [0, 1, 2, 3]:
        x = xpos(float(tick))
        parts.append(f'<line x1="{x}" y1="{top-15}" x2="{x}" y2="{bottom}" stroke="#d8dde5" stroke-width="2"/>')
        parts.append(f'<text x="{x-7}" y="{bottom+38}" class="label">{tick}</text>')
    tx = xpos(THRESHOLD)
    parts.append(f'<line x1="{tx}" y1="{top-25}" x2="{tx}" y2="{bottom}" stroke="#c44e52" stroke-width="6"/>')
    parts.append(f'<text x="{tx-80}" y="{top-45}" class="label bold" style="fill:#c44e52">threshold = 3</text>')

    colors = ["#4c78a8", "#4c78a8", "#7a5dc7", "#7a5dc7"]
    for index, row in enumerate(ROWS):
        y = top + index * 118
        bar_width = xpos(row["ratio"]) - left
        parts.extend(
            [
                f'<text x="60" y="{y+28}" class="label bold">{row["geometry"]}</text>',
                f'<text x="60" y="{y+61}" class="small" style="fill:#596273">contrast {row["contrast"]:.7f} | floor {row["floor"]:.7f}</text>',
                f'<rect x="{left}" y="{y}" width="{bar_width}" height="64" rx="14" fill="{colors[index]}"/>',
                f'<text x="{xpos(row["ratio"])+18}" y="{y+41}" class="label bold">{row["ratio"]:.3f}</text>',
                f'<text x="{right-125}" y="{y+39}" class="small bold" style="fill:#a13f48">UNRESOLVED</text>',
            ]
        )

    parts.extend(
        [
            '<line x1="60" y1="810" x2="1640" y2="810" stroke="#aeb6c2" stroke-width="3"/>',
            '<text x="60" y="865" class="sub">Pooled checks</text>',
            '<text x="60" y="915" class="small">Pooled archived floor</text>',
            '<polygon points="420,900 433,887 446,900 433,913" fill="#4c78a8"/>',
            f'<text x="462" y="907" class="label bold">ratio {POOLED[0][1]:.3f} &lt; 3</text>',
            '<text x="60" y="970" class="small">Pooled conservative floor</text>',
            '<polygon points="950,955 963,942 976,955 963,968" fill="#7a5dc7"/>',
            f'<text x="992" y="962" class="label bold">ratio {POOLED[1][1]:.3f} &lt; 3</text>',
            '<text x="1080" y="905" class="small bold" style="fill:#596273">Below threshold means insufficient audited evidence,</text>',
            '<text x="1080" y="935" class="small bold" style="fill:#596273">not a zero physical interaction.</text>',
            '</svg>',
        ]
    )
    SVG_PATH.write_text("\n".join(parts), encoding="utf-8")


def main() -> None:
    build_png()
    build_svg()
    print(PNG_PATH)
    print(SVG_PATH)


if __name__ == "__main__":
    main()
