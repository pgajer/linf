#!/usr/bin/env python3
"""Build the cross-application synthesis figure for the dCST methods paper."""

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch


ROOT = Path("/Users/pgajer/current_projects/linf/papers/dcst-methods")
OUT_DIR = ROOT / "assets" / "figures"
PDF_OUT = OUT_DIR / "FIGURE_4_cross_application_synthesis.pdf"
PNG_OUT = OUT_DIR / "FIGURE_4_cross_application_synthesis.png"


COLORS = {
    "ink": "#25323a",
    "muted": "#53636d",
    "blue": "#4C78A8",
    "green": "#4F8A55",
    "orange": "#C97A2B",
    "purple": "#7D62A8",
    "paper": "#FBFAF6",
}


def wrapped(text: str, width: int) -> str:
    return "\n".join(textwrap.wrap(text, width=width, break_long_words=False))


def add_header(ax, x, y, w, h, title, subtitle, color):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.5,
        facecolor="#FFFFFF",
        edgecolor=color,
    )
    ax.add_patch(patch)
    ax.text(x + 0.018, y + h - 0.035, title, ha="left", va="top",
            fontsize=12.6, fontweight="bold", color=COLORS["ink"])
    ax.text(x + 0.018, y + h - 0.085, wrapped(subtitle, 36), ha="left", va="top",
            fontsize=8.8, color=COLORS["muted"], linespacing=1.10)


def add_cell(ax, x, y, w, h, title, body, color, wrap=38):
    ax.text(x, y + h - 0.012, title, ha="left", va="top",
            fontsize=11.0, fontweight="bold", color=color)
    ax.text(x, y + h - 0.055, wrapped(body, wrap), ha="left", va="top",
            fontsize=10.0, color=COLORS["ink"], linespacing=1.16)


def add_card(ax, x, y, w, h, title, body, color, wrap=30):
    patch = FancyBboxPatch(
        (x, y),
        w,
        h,
        boxstyle="round,pad=0.010,rounding_size=0.012",
        linewidth=1.0,
        facecolor="#FFFFFF",
        edgecolor="#D7DDD9",
    )
    ax.add_patch(patch)
    ax.text(x + 0.016, y + h - 0.026, title, ha="left", va="top",
            fontsize=10.8, fontweight="bold", color=color)
    ax.text(x + 0.016, y + h - 0.070, wrapped(body, wrap), ha="left", va="top",
            fontsize=9.2, color=COLORS["ink"], linespacing=1.10)


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(8.4, 7.0))
    fig.patch.set_facecolor(COLORS["paper"])
    ax.set_facecolor(COLORS["paper"])
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    ax.text(0.04, 0.965, "Cross-application evaluation of dCSTs",
            ha="left", va="top", fontsize=17.5, fontweight="bold", color=COLORS["ink"])
    ax.text(0.04, 0.915,
            "Two applications, one deterministic dominance-lineage construction.",
            ha="left", va="top", fontsize=10.2, color=COLORS["muted"])

    add_header(
        ax,
        0.04,
        0.755,
        0.43,
        0.120,
        "Gut IBD portability",
        "AGP stool hierarchy tested in external IBD cohorts.",
        COLORS["green"],
    )
    add_header(
        ax,
        0.53,
        0.755,
        0.43,
        0.120,
        "CT full-VOG interpretation",
        "Full-VOG hierarchy tested in labeled CT metagenomes.",
        COLORS["orange"],
    )

    gut_cards = [
        ("State vocabulary", "340 shared stool taxa; AGP frozen hierarchy: 50/158/219/221 states."),
        ("Primary evidence", "Halfvarson retained sample-level IBD, Crohn, and UC signal in transfer and rebuilt modes."),
        ("Boundary", "Subject-aware sensitivity, HMP2, and Gevers bound the portability claim."),
        ("Lesson", "Frozen transfer tests whether a learned taxonomic state vocabulary travels."),
    ]
    ct_cards = [
        ("State vocabulary", "509,501 VOGs; full hierarchy: 12/23/52/68 states."),
        ("Primary evidence", "Depth-4 CT panel: q=0.0047; Cramer's V=0.371 across 58 occupied lineages."),
        ("Boundary", "Feature caps support the panel from 50k VOGs; no independent CT replication."),
        ("Lesson", "Full-VOG dCSTs retain taxonomic, product-level, and KEGG-linked interpretation."),
    ]

    y_positions = [0.590, 0.435, 0.280, 0.125]
    for (title, body), y in zip(gut_cards, y_positions):
        add_card(ax, 0.04, y, 0.43, 0.130, title, body, COLORS["green"], wrap=35)
    for (title, body), y in zip(ct_cards, y_positions):
        add_card(ax, 0.53, y, 0.43, 0.130, title, body, COLORS["orange"], wrap=35)

    ax.text(
        0.04,
        0.050,
        wrapped("Conclusion: dCSTs provide deterministic state vocabularies; portability, robustness, and biology remain empirical properties of each application.", 78),
        ha="left",
        va="bottom",
        fontsize=10.4,
        color=COLORS["purple"],
        fontweight="bold",
    )

    fig.savefig(PNG_OUT, dpi=260, bbox_inches="tight", pad_inches=0.05)
    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.05)


if __name__ == "__main__":
    main()
