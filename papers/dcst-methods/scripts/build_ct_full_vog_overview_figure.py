#!/usr/bin/env python3
"""Build the full-VOG CT overview figure used in the dCST methods paper."""

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path("/Users/pgajer/current_projects/linf/papers/dcst-methods")
OUT_DIR = ROOT / "assets" / "figures"
PDF_OUT = OUT_DIR / "FIGURE_3_ct_full_vog_overview.pdf"
PNG_OUT = OUT_DIR / "FIGURE_3_ct_full_vog_overview.png"


COLORS = {
    "paper": "#FBFAF5",
    "ink": "#263238",
    "muted": "#607D8B",
    "clearance": "#2F80A7",
    "persistence": "#C46A2D",
    "bridge": "#6F4B8B",
    "blue_bg": "#EAF4F7",
    "orange_bg": "#FFF3E6",
    "purple_bg": "#F1E8F8",
}


def wrap(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def box(ax, x, y, w, h, title, body, fc, ec, body_width=40):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.016",
        linewidth=1.25,
        facecolor=fc,
        edgecolor=ec,
    )
    ax.add_patch(patch)
    ax.text(x + 0.018, y + h - 0.022, title, ha="left", va="top",
            fontsize=11.6, fontweight="bold", color=COLORS["ink"])
    ax.text(x + 0.018, y + h - 0.078, wrap(body, body_width), ha="left", va="top",
            fontsize=10.0, color=COLORS["ink"], linespacing=1.18)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 6.8))
    fig.patch.set_facecolor(COLORS["paper"])
    ax.set_facecolor(COLORS["paper"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.035, 0.955, "Full 509,501-VOG dCST-CT association layer",
            fontsize=18, fontweight="bold", color=COLORS["ink"], va="top")
    ax.text(0.035, 0.895,
            wrap("Outcome tests use 604 labeled CT metagenomes; the hierarchy is learned in the 3,602-sample CT+reference cohort.", 70),
            fontsize=10.2, color=COLORS["muted"], va="top", linespacing=1.08)

    box(
        ax,
        0.035,
        0.670,
        0.93,
        0.145,
        "Panel-level result",
        "Depth-4 association: omnibus q=0.0047; Cramer's V=0.371. Capped runs retain the panel signal from 50k VOGs onward.",
        "#E6F0FA",
        "#315F8B",
        body_width=56,
    )

    rows = [
        (
            0.445,
            "Clearance branch",
            "18 CT: 15 clearance / 3 persistence. Bayesian OR=0.20 (0.05-0.57).",
            "Mostly L. jensenii; L. crispatus root.",
            "Sugar/ABC transport; D-lactate dehydrogenase; cation ATPases.",
            COLORS["blue_bg"],
            COLORS["clearance"],
        ),
        (
            0.250,
            "Persistence branch",
            "20 CT: 4 clearance / 16 persistence. Bayesian OR=3.85 (1.38-12.78).",
            "Lactobacillus iners.",
            "ABC transport; threonyl-tRNA ligase; regulators; surface-anchor proteins.",
            COLORS["orange_bg"],
            COLORS["persistence"],
        ),
        (
            0.055,
            "Persistence bridge",
            "31 CT: 9 clearance / 22 persistence. Bayesian OR=2.46 (1.15-5.72).",
            "Ca. Lachnocurva vaginae.",
            "Tungsten/ABC transport; phosphatase; motility-linked proteins.",
            COLORS["purple_bg"],
            COLORS["bridge"],
        ),
    ]

    for y, signal_title, signal, taxonomy, function, fc, ec in rows:
        box(ax, 0.035, y, 0.305, 0.165, signal_title, signal, fc, ec, body_width=25)
        box(ax, 0.360, y, 0.245, 0.165, "Taxonomy", taxonomy, "#FFFFFF", ec, body_width=18)
        box(ax, 0.625, y, 0.340, 0.165, "Function", function, "#FFFFFF", ec, body_width=27)

    fig.savefig(PNG_OUT, dpi=260, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.05)


if __name__ == "__main__":
    main()
