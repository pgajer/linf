from __future__ import annotations

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import FancyArrowPatch, Rectangle

from paper_paths import FIGURES_DIR
from taxon_formatting import italicize_taxa_mpl

OUT = FIGURES_DIR / "FIGURE_1_dcst_conceptual_schematic.png"


def add_panel_label(ax, label: str) -> None:
    ax.text(
        -0.08,
        1.05,
        label,
        transform=ax.transAxes,
        fontsize=13,
        fontweight="bold",
        va="top",
        ha="left",
    )


def draw_box(ax, xy, w, h, text, face, edge="#4f4f4f", fontsize=8.0, weight=None):
    patch = Rectangle(xy, w, h, facecolor=face, edgecolor=edge, linewidth=1.2)
    ax.add_patch(patch)
    ax.text(
        xy[0] + w / 2,
        xy[1] + h / 2,
        text,
        ha="center",
        va="center",
        fontsize=fontsize,
        fontweight=weight,
        linespacing=1.05,
    )
    return patch


def arrow(ax, start, end):
    ax.add_patch(
        FancyArrowPatch(
            start,
            end,
            arrowstyle="-|>",
            mutation_scale=14,
            color="#6c757d",
            linewidth=1.2,
        )
    )


def italic_label(label: str) -> str:
    return italicize_taxa_mpl(label)


def italic_two_line(first: str, second: str) -> str:
    return italic_label(first) + "-\n" + italic_label(second)


def main() -> None:
    plt.rcParams.update(
        {
            "font.size": 8.0,
            "axes.titlesize": 8.4,
            "axes.labelsize": 7.5,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    taxa = ["Bacteroides", "Faecalibacterium", "Escherichia-Shigella", "Alistipes"]
    samples = [f"S{i}" for i in range(1, 7)]
    matrix = np.array(
        [
            [0.48, 0.22, 0.10, 0.20],
            [0.42, 0.31, 0.09, 0.18],
            [0.38, 0.13, 0.16, 0.33],
            [0.18, 0.49, 0.12, 0.21],
            [0.15, 0.40, 0.28, 0.17],
            [0.20, 0.12, 0.48, 0.20],
        ]
    )

    fig = plt.figure(figsize=(7.0, 5.25))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.0], hspace=0.72, wspace=0.36)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, 0])
    ax_d = fig.add_subplot(gs[1, 1])
    fig.subplots_adjust(left=0.075, right=0.985, top=0.94, bottom=0.125, hspace=0.74, wspace=0.34)

    im = ax_a.imshow(matrix, aspect="auto", cmap="YlOrRd", vmin=0, vmax=0.5)
    ax_a.set_xticks(np.arange(len(taxa)))
    ax_a.set_xticklabels([italic_label(t) for t in taxa], rotation=28, ha="right", fontsize=6.4)
    ax_a.set_yticks(np.arange(len(samples)))
    ax_a.set_yticklabels(samples, fontsize=6.6)
    ax_a.set_title("Toy gut abundance table", loc="left", fontweight="bold")
    ax_a.set_xlabel("Taxa")
    ax_a.set_ylabel("Samples")
    add_panel_label(ax_a, "A")
    cbar = fig.colorbar(im, ax=ax_a, fraction=0.046, pad=0.02)
    cbar.set_label("Relative abundance", fontsize=7.0)

    ax_b.axis("off")
    ax_b.set_xlim(0, 1)
    ax_b.set_ylim(0, 1)
    add_panel_label(ax_b, "B")
    ax_b.set_title("Depth-1 dominance sample sets", loc="left", fontweight="bold")
    draw_box(ax_b, (0.01, 0.69), 0.39, 0.18, "Top taxon\n" + italic_label("Bacteroides"), "#d9ead3", "#6a994e")
    draw_box(ax_b, (0.01, 0.44), 0.39, 0.18, "Top taxon\n" + italic_label("Faecalibacterium"), "#cfe2f3", "#457b9d")
    draw_box(ax_b, (0.01, 0.17), 0.39, 0.20, "Top taxon\n" + italic_two_line("Escherichia", "Shigella"), "#fce5cd", "#bc6c25", fontsize=7.2)
    draw_box(ax_b, (0.47, 0.69), 0.52, 0.18, "Bacteroides DSS\nn = 3", "#d9ead3", "#6a994e", weight="bold")
    draw_box(ax_b, (0.47, 0.44), 0.52, 0.18, "Faecalibacterium DSS\nn = 2", "#cfe2f3", "#457b9d", weight="bold")
    draw_box(ax_b, (0.47, 0.17), 0.52, 0.20, "Escherichia-\nShigella DSS\nn = 1", "#fce5cd", "#bc6c25", fontsize=6.8, weight="bold")
    for y in (0.78, 0.53, 0.27):
        arrow(ax_b, (0.40, y), (0.47, y))

    ax_c.axis("off")
    ax_c.set_xlim(0, 1)
    ax_c.set_ylim(0, 1)
    add_panel_label(ax_c, "C")
    ax_c.set_title("Low-frequency policy", loc="left", fontweight="bold")
    draw_box(ax_c, (0.03, 0.62), 0.36, 0.20, "Small DSS\nn < n0", "#fff2cc", "#bf9000")
    draw_box(ax_c, (0.47, 0.73), 0.51, 0.17, "Pure dCST:\nrare/residual bucket", "#eeeeee", "#6c757d", fontsize=7.4, weight="bold")
    draw_box(ax_c, (0.47, 0.43), 0.51, 0.18, "Absorb dCST:\nnearest retained state", "#eaf4f4", "#4d908e", fontsize=7.1, weight="bold")
    arrow(ax_c, (0.39, 0.74), (0.47, 0.815))
    arrow(ax_c, (0.39, 0.67), (0.47, 0.52))
    draw_box(ax_c, (0.04, 0.19), 0.40, 0.16, "Retained DSS\nn >= n0", "#e2f0d9", "#6a994e")
    draw_box(ax_c, (0.52, 0.19), 0.41, 0.16, "Named dCST", "#edf6f9", "#457b9d", fontsize=7.8, weight="bold")
    arrow(ax_c, (0.44, 0.27), (0.52, 0.27))

    ax_d.axis("off")
    ax_d.set_xlim(0, 1)
    ax_d.set_ylim(0, 1)
    add_panel_label(ax_d, "D")
    ax_d.set_title("Depth-2 dCST dominance-lineages", loc="left", fontweight="bold")
    draw_box(ax_d, (0.03, 0.67), 0.43, 0.18, "Parent DSS\n" + italic_label("Bacteroides"), "#d9ead3", "#6a994e", weight="bold")
    draw_box(ax_d, (0.52, 0.76), 0.47, 0.16, italic_label("Bacteroides") + " /\n" + italic_label("Faecalibacterium"), "#cfe2f3", "#457b9d", fontsize=7.5)
    draw_box(ax_d, (0.52, 0.53), 0.47, 0.16, italic_label("Bacteroides") + " /\n" + italic_label("Alistipes"), "#ead1dc", "#8e7cc3", fontsize=7.5)
    arrow(ax_d, (0.46, 0.76), (0.52, 0.84))
    arrow(ax_d, (0.46, 0.74), (0.52, 0.61))
    draw_box(ax_d, (0.06, 0.18), 0.91, 0.20, "Ordered path:\nparent dominant taxon / next dominant taxon", "#f6f1e8", "#b08968", fontsize=7.0)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=340, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    main()
