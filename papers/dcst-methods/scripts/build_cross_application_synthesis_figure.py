#!/usr/bin/env python3

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
from PIL import Image


ROOT = Path("/Users/pgajer/current_projects/linf/papers/dcst-methods")
OUT_DIR = ROOT / "assets" / "figures"
PDF_OUT = OUT_DIR / "FIGURE_5_cross_application_synthesis.pdf"
PNG_OUT = OUT_DIR / "FIGURE_5_cross_application_synthesis.png"


def add_box(
    ax,
    xy,
    width,
    height,
    facecolor,
    edgecolor,
    title,
    lines,
    title_color="#24303a",
    title_size=15,
    body_size=12.2,
    body_y_offset=0.105,
    wrap_width=None,
):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.02",
        linewidth=1.8,
        facecolor=facecolor,
        edgecolor=edgecolor,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.02,
        y + height - 0.055,
        title,
        ha="left",
        va="top",
        fontsize=title_size,
        fontweight="bold",
        color=title_color,
    )
    if wrap_width is None:
        body = "\n".join(lines)
    else:
        body = "\n".join(
            textwrap.fill(line, width=wrap_width, break_long_words=False)
            for line in lines
        )
    ax.text(
        x + 0.02,
        y + height - body_y_offset,
        body,
        ha="left",
        va="top",
        fontsize=body_size,
        color="#24303a",
        linespacing=1.35,
    )


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10.6, 7.5))
    fig.patch.set_facecolor("white")
    ax.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    add_box(
        ax,
        (0.06, 0.735),
        0.88,
        0.205,
        facecolor="#eef4fb",
        edgecolor="#4c78a8",
        title="Shared dCST framework",
        lines=[
            "Retained feature family + depth-wise dominance-lineage recursion",
            "Absorb-policy labels separate state construction from downstream inference",
        ],
        title_size=16,
        body_size=10.8,
        body_y_offset=0.098,
        wrap_width=72,
    )

    add_box(
        ax,
        (0.06, 0.42),
        0.41,
        0.26,
        facecolor="#eef8f1",
        edgecolor="#4f8a55",
        title="Gut IBD portability",
        lines=[
            "AGP-derived shared-taxonomy states transferred and rebuilt in external IBD cohorts",
            "Signal is clinically coherent in Halfvarson but attenuates under subject-aware sensitivity",
        ],
        title_size=16,
        body_size=11.8,
        wrap_width=36,
    )

    add_box(
        ax,
        (0.53, 0.42),
        0.41,
        0.26,
        facecolor="#fbf1e8",
        edgecolor="#c97a2b",
        title="CT functional interpretation",
        lines=[
            "VOG-cluster and full-VOG states associate with clearance and persistence",
            "Capped VOG robustness and resolution-level concordance support panel-level organization",
        ],
        title_size=16,
        body_size=11.8,
        wrap_width=36,
    )

    add_box(
        ax,
        (0.06, 0.08),
        0.88,
        0.22,
        facecolor="#f7f4fa",
        edgecolor="#7d62a8",
        title="Cross-application synthesis",
        lines=[
            "Deterministic dominance-lineages provide a reusable microbiome state representation.",
            "Portability, robustness, and biological interpretability remain application-dependent empirical properties.",
        ],
        title_size=16,
        body_size=11.8,
        body_y_offset=0.083,
        wrap_width=70,
    )

    arrow_style = dict(arrowstyle="-|>", linewidth=1.8, color="#6b7b8c", mutation_scale=14)
    ax.annotate("", xy=(0.265, 0.68), xytext=(0.265, 0.735), arrowprops=arrow_style)
    ax.annotate("", xy=(0.735, 0.68), xytext=(0.735, 0.735), arrowprops=arrow_style)
    ax.annotate("", xy=(0.50, 0.30), xytext=(0.50, 0.42), arrowprops=arrow_style)

    fig.tight_layout()
    fig.savefig(PNG_OUT, dpi=220, bbox_inches="tight")
    with Image.open(PNG_OUT) as img:
        img.convert("RGB").save(PDF_OUT, "PDF", resolution=220.0)


if __name__ == "__main__":
    main()
