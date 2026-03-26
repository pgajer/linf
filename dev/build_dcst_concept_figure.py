from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle


OUT = Path("/Users/pgajer/current_projects/linf/dev/FIGURE_1_dcst_conceptual_schematic.png")


def draw_bars(ax, values, colors, title, xlabel_text, annotate=None):
    x = np.arange(len(values))
    ax.bar(x, values, color=colors, edgecolor="#2f2f2f", linewidth=0.8)
    ax.set_xticks(x)
    ax.set_xticklabels(["Taxon A", "Taxon B", "Taxon C", "Taxon D"], rotation=20, ha="right")
    ax.set_ylim(0, max(1.08, max(values) * 1.18))
    ax.set_title(title, fontsize=12, pad=10, weight="bold")
    ax.set_xlabel(xlabel_text, fontsize=10)
    ax.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    if annotate:
        for idx, txt in annotate:
            ax.text(idx, values[idx] + 0.04, txt, ha="center", va="bottom", fontsize=9)


def add_panel_label(ax, label):
    ax.text(
        -0.13,
        1.08,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def main():
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    raw = np.array([420, 280, 180, 120])
    colors = ["#335c67", "#7cb518", "#ff9f1c", "#c8553d"]

    fig = plt.figure(figsize=(13, 8))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.08], hspace=0.42, wspace=0.28)

    ax1 = fig.add_subplot(gs[0, 0])
    draw_bars(
        ax1,
        raw,
        colors,
        "Sample abundance profile",
        r"Example sample vector $x = (420,\, 280,\, 180,\, 120)$",
    )
    add_panel_label(ax1, "A")

    ax2 = fig.add_subplot(gs[0, 1])
    draw_bars(
        ax2,
        raw,
        colors,
        "Within-sample dominance order",
        "Only the rank order matters: Taxon A > Taxon B > Taxon C > Taxon D",
        annotate=[(0, "1st"), (1, "2nd"), (2, "3rd"), (3, "4th")],
    )
    add_panel_label(ax2, "B")

    ax3 = fig.add_subplot(gs[1, 0])
    ax3.axis("off")
    add_panel_label(ax3, "C")
    box_fc = "#f6f1e8"
    box_ec = "#b08968"
    text_color = "#2f2f2f"

    left = Rectangle((0.03, 0.52), 0.40, 0.26, facecolor=box_fc, edgecolor=box_ec, linewidth=1.5)
    right = Rectangle((0.56, 0.52), 0.37, 0.26, facecolor="#edf6f9", edgecolor="#4d908e", linewidth=1.5)
    ax3.add_patch(left)
    ax3.add_patch(right)
    ax3.text(0.23, 0.65, "Largest coordinate\nTaxon A", ha="center", va="center", fontsize=12, color=text_color)
    ax3.text(0.745, 0.65, "Depth-1 DCST\nA", ha="center", va="center", fontsize=13, weight="bold", color="#1d3557")
    arrow = FancyArrowPatch((0.43, 0.65), (0.56, 0.65), arrowstyle="simple", mutation_scale=18, color="#6c757d")
    ax3.add_patch(arrow)
    ax3.text(
        0.5,
        0.25,
        "Depth 1 keeps only the top-ranked taxon.\nThe label is read directly from the sample, not from a fitted cluster.",
        ha="center",
        va="center",
        fontsize=10.5,
        color=text_color,
    )

    ax4 = fig.add_subplot(gs[1, 1])
    ax4.axis("off")
    add_panel_label(ax4, "D")
    box1 = Rectangle((0.03, 0.57), 0.25, 0.22, facecolor=box_fc, edgecolor=box_ec, linewidth=1.5)
    box2 = Rectangle((0.32, 0.57), 0.25, 0.22, facecolor="#e9f5db", edgecolor="#6a994e", linewidth=1.5)
    box3 = Rectangle((0.62, 0.57), 0.32, 0.22, facecolor="#edf6f9", edgecolor="#4d908e", linewidth=1.5)
    for patch in (box1, box2, box3):
        ax4.add_patch(patch)
    ax4.text(0.155, 0.68, "Top rank\nTaxon A", ha="center", va="center", fontsize=11.5, color=text_color)
    ax4.text(0.445, 0.68, "Second rank\nTaxon B", ha="center", va="center", fontsize=11.5, color=text_color)
    ax4.text(0.78, 0.68, "Depth-2 DCST\nA__B", ha="center", va="center", fontsize=13, weight="bold", color="#1d3557")
    ax4.add_patch(FancyArrowPatch((0.28, 0.68), (0.32, 0.68), arrowstyle="-|>", mutation_scale=16, color="#6c757d"))
    ax4.add_patch(FancyArrowPatch((0.57, 0.68), (0.62, 0.68), arrowstyle="-|>", mutation_scale=16, color="#6c757d"))
    ax4.text(
        0.48,
        0.27,
        "Depth 2 records the ordered dominant pair.\nIn the rare policy used here, low-frequency patterns are grouped into RARE_DOMINANT\nrather than absorbed into larger named states.",
        ha="center",
        va="center",
        fontsize=10.5,
        color=text_color,
    )

    fig.suptitle(
        "Figure 1. Conceptual construction of dominant-community state types (DCSTs)",
        fontsize=15,
        weight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.02,
        "The key point is that the label is determined by within-sample dominance order, yielding a deterministic and hierarchical alternative to clustering-based community typing.",
        ha="center",
        va="bottom",
        fontsize=10.5,
        color="#3b3b3b",
    )
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=220, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
