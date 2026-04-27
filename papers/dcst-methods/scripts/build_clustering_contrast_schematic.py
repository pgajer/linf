#!/usr/bin/env python3

from pathlib import Path
import textwrap

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from PIL import Image


ROOT = Path("/Users/pgajer/current_projects/linf/papers/dcst-methods")
OUT_DIR = ROOT / "assets" / "figures"
PNG_OUT = OUT_DIR / "FIGURE_S2_clustering_contrast_schematic.png"
PDF_OUT = OUT_DIR / "FIGURE_S2_clustering_contrast_schematic.pdf"


def add_box(ax, xy, width, height, title, body, face, edge, body_size=10.5):
    x, y = xy
    patch = FancyBboxPatch(
        (x, y),
        width,
        height,
        boxstyle="round,pad=0.012,rounding_size=0.018",
        linewidth=1.5,
        facecolor=face,
        edgecolor=edge,
    )
    ax.add_patch(patch)
    ax.text(
        x + 0.018,
        y + height - 0.045,
        title,
        ha="left",
        va="top",
        fontsize=13,
        fontweight="bold",
        color="#24303a",
    )
    ax.text(
        x + 0.018,
        y + height - 0.088,
        "\n".join(textwrap.wrap(body, width=58, break_long_words=False)),
        ha="left",
        va="top",
        fontsize=body_size,
        color="#24303a",
        linespacing=1.25,
    )


def add_label_row(ax, y, title, labels, colors, x0=0.08, dx=0.055, title_x=None, label_width=0.043):
    if title_x is None:
        title_x = x0 - 0.035
    ax.text(title_x, y, title, ha="right", va="center", fontsize=8.8, color="#24303a", linespacing=1.0)
    for i, (label, color) in enumerate(zip(labels, colors)):
        x = x0 + i * dx
        ax.add_patch(
            FancyBboxPatch(
                (x, y - 0.017),
                label_width,
                0.034,
                boxstyle="round,pad=0.006,rounding_size=0.008",
                linewidth=0.8,
                facecolor=color,
                edgecolor="#39424c",
            )
        )
        ax.text(x + label_width / 2, y, label, ha="center", va="center", fontsize=8.0, color="#111820")


def add_sample_ids(ax, y, ids, x0=0.08, dx=0.055, title_x=None, label_width=0.043):
    if title_x is None:
        title_x = x0 - 0.035
    ax.text(title_x, y, "Target\nsamples", ha="right", va="center", fontsize=8.4, color="#56616f", linespacing=1.0)
    for i, sample_id in enumerate(ids):
        ax.text(x0 + i * dx + label_width / 2, y, sample_id, ha="center", va="center", fontsize=8.1, color="#56616f")


def add_dots(ax, xs, ys, colors, labels=None):
    for i, (x, y, color) in enumerate(zip(xs, ys, colors)):
        ax.scatter([x], [y], s=78, color=color, edgecolor="#25313d", linewidth=0.6, zorder=3)
        if labels:
            ax.text(x, y - 0.038, labels[i], ha="center", va="top", fontsize=8.0, color="#47515c")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    fig, ax = plt.subplots(figsize=(10.8, 6.7))
    fig.patch.set_facecolor("white")
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.axis("off")

    colors = {
        "blue": "#8fb9dd",
        "green": "#9fd0a8",
        "orange": "#f1b273",
        "purple": "#c9b4e3",
        "gray": "#d9dde2",
    }

    ax.text(
        0.04,
        0.955,
        "Toy schematic: clustering labels are conditional; frozen dCST assignments are fixed",
        ha="left",
        va="top",
        fontsize=14.2,
        fontweight="bold",
        color="#24303a",
    )
    ax.text(
        0.04,
        0.915,
        "This schematic illustrates label definition only; it is not a benchmark of clustering performance.",
        ha="left",
        va="top",
        fontsize=9.8,
        color="#56616f",
    )

    left = Rectangle((0.035, 0.12), 0.45, 0.75, facecolor="#f7f8fa", edgecolor="#c9d1db", linewidth=1.1)
    right = Rectangle((0.515, 0.12), 0.45, 0.75, facecolor="#f8faf7", edgecolor="#cbd8ca", linewidth=1.1)
    ax.add_patch(left)
    ax.add_patch(right)

    ax.text(0.055, 0.835, "A. Clustering-defined CST labels", ha="left", va="top", fontsize=12.4, fontweight="bold", color="#24303a")
    ax.text(0.535, 0.835, "B. Frozen dCST assignment", ha="left", va="top", fontsize=12.4, fontweight="bold", color="#24303a")

    ids = ["s1", "s2", "s3", "s4", "s5", "s6"]
    add_sample_ids(ax, 0.765, ids, x0=0.15, dx=0.052, title_x=0.115)
    add_label_row(
        ax,
        0.705,
        "Cluster,\nk = 2",
        ["C1", "C1", "C1", "C2", "C2", "C2"],
        [colors["blue"]] * 3 + [colors["orange"]] * 3,
        x0=0.15,
        dx=0.052,
        title_x=0.115,
    )
    add_label_row(
        ax,
        0.645,
        "Cluster,\nk = 3",
        ["C1", "C1", "C2", "C2", "C3", "C3"],
        [colors["blue"]] * 2 + [colors["green"]] * 2 + [colors["orange"]] * 2,
        x0=0.15,
        dx=0.052,
        title_x=0.115,
    )
    add_label_row(
        ax,
        0.585,
        "Add\nsamples",
        ["C1", "C2", "C2", "C3", "C3", "C3"],
        [colors["blue"]] + [colors["green"]] * 2 + [colors["orange"]] * 3,
        x0=0.15,
        dx=0.052,
        title_x=0.115,
    )

    ax.plot([0.075, 0.42], [0.535, 0.535], color="#c9d1db", linewidth=1)
    add_box(
        ax,
        (0.075, 0.315),
        0.37,
        0.165,
        "Interpretation",
        "The label is conditional on choices such as k, distance, algorithm, initialization, and the comparison set used to fit the partition.",
        face="#ffffff",
        edge="#c9d1db",
        body_size=8.8,
    )
    ax.text(0.075, 0.245, "Same target samples; partitions can differ.", ha="left", va="top", fontsize=8.9, color="#56616f")

    ax.text(0.545, 0.765, "Frozen reference hierarchy", ha="left", va="center", fontsize=10.5, fontweight="bold", color="#24303a")
    ax.add_patch(FancyBboxPatch((0.56, 0.675), 0.13, 0.055, boxstyle="round,pad=0.008", facecolor="#e9f2fb", edgecolor="#4c78a8"))
    ax.add_patch(FancyBboxPatch((0.75, 0.675), 0.13, 0.055, boxstyle="round,pad=0.008", facecolor="#fdf1e5", edgecolor="#c97a2b"))
    ax.text(0.625, 0.702, "A parent", ha="center", va="center", fontsize=9.5, fontweight="bold", color="#24303a")
    ax.text(0.815, 0.702, "B parent", ha="center", va="center", fontsize=9.5, fontweight="bold", color="#24303a")
    for x, label, color in [(0.535, "A>B", colors["blue"]), (0.625, "A>C", colors["green"]), (0.715, "B>A", colors["orange"]), (0.805, "B>C", colors["purple"])]:
        ax.add_patch(FancyBboxPatch((x, 0.595), 0.07, 0.045, boxstyle="round,pad=0.006", facecolor=color, edgecolor="#39424c", linewidth=0.8))
        ax.text(x + 0.035, 0.6175, label, ha="center", va="center", fontsize=9.2, color="#111820")
    ax.annotate("", xy=(0.57, 0.642), xytext=(0.61, 0.675), arrowprops=dict(arrowstyle="-", linewidth=1.0, color="#607080"))
    ax.annotate("", xy=(0.66, 0.642), xytext=(0.64, 0.675), arrowprops=dict(arrowstyle="-", linewidth=1.0, color="#607080"))
    ax.annotate("", xy=(0.75, 0.642), xytext=(0.79, 0.675), arrowprops=dict(arrowstyle="-", linewidth=1.0, color="#607080"))
    ax.annotate("", xy=(0.84, 0.642), xytext=(0.83, 0.675), arrowprops=dict(arrowstyle="-", linewidth=1.0, color="#607080"))

    add_sample_ids(ax, 0.51, ids, x0=0.64, dx=0.047, title_x=0.61, label_width=0.039)
    add_label_row(
        ax,
        0.45,
        "Before\nadding",
        ["A>B", "A>B", "A>C", "B>A", "B>C", "A>C"],
        [colors["blue"], colors["blue"], colors["green"], colors["orange"], colors["purple"], colors["green"]],
        x0=0.64,
        dx=0.047,
        title_x=0.61,
        label_width=0.039,
    )
    add_label_row(
        ax,
        0.39,
        "After\nadding",
        ["A>B", "A>B", "A>C", "B>A", "B>C", "A>C"],
        [colors["blue"], colors["blue"], colors["green"], colors["orange"], colors["purple"], colors["green"]],
        x0=0.64,
        dx=0.047,
        title_x=0.61,
        label_width=0.039,
    )
    add_box(
        ax,
        (0.545, 0.195),
        0.37,
        0.125,
        "Interpretation",
        "Frozen hierarchy: aligned samples map to the same realized lineages.",
        face="#ffffff",
        edge="#cbd8ca",
        body_size=8.8,
    )

    fig.savefig(PNG_OUT, dpi=240, bbox_inches="tight")
    with Image.open(PNG_OUT) as img:
        img.convert("RGB").save(PDF_OUT, "PDF", resolution=240.0)


if __name__ == "__main__":
    main()
