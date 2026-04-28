#!/usr/bin/env python3
"""Build an alternate Figure 1 with manuscript-consistent dCST labels.

This script intentionally writes a new asset and does not overwrite the
original FIGURE_1_dcst_conceptual_framework.pdf.
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import FancyArrowPatch


ROOT = Path("/Users/pgajer/current_projects/linf/papers/dcst-methods")
OUT_DIR = ROOT / "assets" / "figures"
PDF_OUT = OUT_DIR / "FIGURE_1_dcst_conceptual_framework_text_labels.pdf"
PNG_OUT = OUT_DIR / "FIGURE_1_dcst_conceptual_framework_text_labels.png"

CMAP = LinearSegmentedColormap.from_list("dcst_yellow_red", ["#fff600", "#ff9a00", "#e60000"])
INK = "#1f2933"
MUTED = "#4a5560"
ARROW = "#4a4a4a"


def make_block(rows, cols, dominant_cols, noise=0.08, seed=1):
    rng = np.random.default_rng(seed)
    mat = rng.uniform(0.0, noise, size=(rows, cols))
    for col, strength in dominant_cols:
        mat[:, col] += strength
    return np.clip(mat, 0, 1)


def stack_blocks(blocks, gap=1):
    width = max(block.shape[1] for block in blocks)
    padded = []
    for block in blocks:
        if block.shape[1] < width:
            pad = np.full((block.shape[0], width - block.shape[1]), np.nan)
            block = np.hstack([block, pad])
        padded.append(block)
        padded.append(np.full((gap, width), np.nan))
    return np.vstack(padded[:-1])


def add_heatmap(fig, rect, data, features, title=None, right_labels=None, label_size=7.5):
    ax = fig.add_axes(rect)
    masked = np.ma.masked_invalid(data)
    ax.imshow(masked, cmap=CMAP, vmin=0, vmax=1, aspect="auto", interpolation="nearest")
    ax.set_xticks(np.arange(len(features)))
    ax.set_xticklabels(features, fontsize=label_size)
    ax.set_yticks([])
    ax.tick_params(axis="x", length=2, pad=1)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if title:
        ax.set_title(title, fontsize=9.0, fontweight="bold", color=INK, pad=7)
    if right_labels:
        total = data.shape[0]
        for label, y_frac in right_labels:
            ax.text(1.02, y_frac, label, transform=ax.transAxes, ha="left", va="center",
                    fontsize=label_size, fontweight="bold", color="black")
    return ax


def add_note(fig, x, y, text, color=INK, size=7.2, weight="normal", ha="left"):
    fig.text(x, y, text, ha=ha, va="top", fontsize=size, color=color,
             fontweight=weight, linespacing=1.12)


def add_arrow(fig, start, end, text=None, rotation=0, size=10):
    arrow = FancyArrowPatch(
        start,
        end,
        transform=fig.transFigure,
        arrowstyle="-|>",
        mutation_scale=10,
        linewidth=1.2,
        color=ARROW,
    )
    fig.patches.append(arrow)
    if text:
        mx = (start[0] + end[0]) / 2
        my = (start[1] + end[1]) / 2
        fig.text(mx, my + 0.012, text, ha="center", va="center",
                 fontsize=size, color=MUTED, rotation=rotation)


def add_panel_label(fig, x, y, label):
    fig.text(x, y, label, ha="left", va="top", fontsize=13,
             fontweight="bold", color="black")


def main():
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    fig = plt.figure(figsize=(9.2, 11.2), facecolor="white")

    add_panel_label(fig, 0.055, 0.972, "A. Depth-1 dCST construction")

    features8 = [f"f{i}" for i in range(1, 9)]
    f1 = make_block(38, 8, [(0, 0.62), (1, 0.30)], seed=10)
    f4 = make_block(32, 8, [(3, 0.70), (4, 0.34), (5, 0.24)], seed=11)
    f7 = make_block(6, 8, [(6, 0.88)], seed=12)
    f8 = make_block(5, 8, [(7, 0.82), (5, 0.35)], seed=13)
    input_a = stack_blocks([f1, f4, f7, f8], gap=1)

    add_heatmap(fig, [0.07, 0.58, 0.20, 0.34], input_a, features8, "Input abundance matrix")
    fig.text(0.17, 0.555, "Feature", ha="center", fontsize=8, color=INK)
    add_arrow(fig, (0.29, 0.765), (0.39, 0.765), "identify dominant\nsample sets", size=7.5)

    add_heatmap(
        fig,
        [0.41, 0.61, 0.16, 0.31],
        input_a,
        features8,
        "Depth-1 dominance\nsample sets",
        [
            ("dominant f1", 0.80),
            ("dominant f4", 0.43),
            ("dominant f7", 0.12),
            ("dominant f8", 0.03),
        ],
        label_size=7,
    )
    fig.text(0.49, 0.585, "Feature", ha="center", fontsize=8, color=INK)

    add_arrow(fig, (0.60, 0.765), (0.74, 0.855), "absorb low-support\nbranches", rotation=36, size=7)
    add_arrow(fig, (0.60, 0.740), (0.74, 0.635), "keep residual", rotation=-36, size=7)

    absorb_a = stack_blocks([np.vstack([f1, f7]), np.vstack([f4, f8])], gap=2)
    add_heatmap(
        fig,
        [0.75, 0.765, 0.15, 0.135],
        absorb_a,
        features8,
        "Absorb policy:\nretained dCST states",
        None,
        label_size=6.5,
    )
    add_note(fig, 0.905, 0.858, "f1 state = dominant f1\n+ absorbed f7", size=6.4, weight="bold")
    add_note(fig, 0.905, 0.782, "f4 state = dominant f4\n+ absorbed f8", size=6.4, weight="bold")
    fig.text(0.84, 0.735, "Feature", ha="center", fontsize=8, color=INK)

    pure_a = stack_blocks([f1, f4, np.vstack([f7, f8])], gap=2)
    add_heatmap(
        fig,
        [0.75, 0.545, 0.15, 0.155],
        pure_a,
        features8,
        "Pure policy:\nretained states + residual",
        None,
        label_size=6.5,
    )
    add_note(fig, 0.905, 0.655, "f1 state", size=6.4, weight="bold")
    add_note(fig, 0.905, 0.607, "f4 state", size=6.4, weight="bold")
    add_note(fig, 0.905, 0.565, "residual = dominant\nf7 or f8", size=6.4, weight="bold")
    fig.text(0.84, 0.520, "Feature", ha="center", fontsize=8, color=INK)

    add_panel_label(fig, 0.055, 0.465, "B. Depth-2 dCST construction")
    fig.text(0.16, 0.425, "Parent dCST state", ha="center", fontsize=8.5, fontweight="bold", color=INK)
    fig.text(0.39, 0.425, "Pre-policy depth-2\nlineages", ha="center", fontsize=8.5, fontweight="bold", color=INK)
    fig.text(0.59, 0.425, "Absorb-policy retained\ndepth-2 states", ha="center", fontsize=8.5, fontweight="bold", color=INK)
    fig.text(0.82, 0.425, "Pure-policy states\nplus residual", ha="center", fontsize=8.5, fontweight="bold", color=INK)

    def depth2_lane(y, parent_name, features, seeds, labels):
        fig.text(0.055, y + 0.140, f"Refine {parent_name} parent", ha="left",
                 fontsize=8.6, fontweight="bold", color="black")
        parent = stack_blocks([
            make_block(16, 4, [(0, 0.58), (1, 0.22)], seed=seeds[0]),
            make_block(5, 4, [(2, 0.75)], seed=seeds[1]),
            make_block(4, 4, [(3, 0.78)], seed=seeds[2]),
        ], gap=1)
        pre1 = make_block(16, 4, [(0, 0.55), (1, 0.34)], seed=seeds[3])
        pre2 = make_block(5, 4, [(2, 0.75)], seed=seeds[4])
        pre3 = make_block(4, 4, [(3, 0.78)], seed=seeds[5])
        pre = stack_blocks([pre1, pre2, pre3], gap=1)
        absorb = stack_blocks([np.vstack([pre1, pre3]), pre2], gap=1)
        pure = pre

        add_heatmap(fig, [0.10, y, 0.13, 0.125], parent, features, None, None, label_size=6.2)
        add_arrow(fig, (0.25, y + 0.065), (0.32, y + 0.065), "refine by next\nhighest feature", size=6.5)
        add_heatmap(fig, [0.35, y, 0.13, 0.125], pre, features, None, None, label_size=6.2)
        add_heatmap(fig, [0.55, y, 0.13, 0.125], absorb, features, None, None, label_size=6.2)
        add_heatmap(fig, [0.78, y, 0.13, 0.125], pure, features, None, None, label_size=6.2)
        fig.text(0.165, y - 0.022, "Feature", ha="center", fontsize=6.5, color=INK)
        fig.text(0.415, y - 0.022, "Feature", ha="center", fontsize=6.5, color=INK)
        fig.text(0.615, y - 0.022, "Feature", ha="center", fontsize=6.5, color=INK)
        fig.text(0.845, y - 0.022, "Feature", ha="center", fontsize=6.5, color=INK)

    depth2_lane(
        0.265,
        "f1",
        ["f1", "f2", "f3", "f9"],
        [20, 21, 22, 23, 24, 25],
        ["f1 > f2 lineage", "f1 > f3 lineage", "f1 > f9 lineage"],
    )
    depth2_lane(
        0.095,
        "f4",
        ["f4", "f5", "f6", "f10"],
        [30, 31, 32, 33, 34, 35],
        ["f4 > f5 lineage", "f4 > f6 lineage", "f4 > f10 lineage"],
    )

    fig.text(0.50, 0.010, "Feature abundance scale: yellow = 0, red = 1.",
             ha="center", fontsize=10, style="italic", color=INK)

    fig.savefig(PDF_OUT, bbox_inches="tight", pad_inches=0.12)
    fig.savefig(PNG_OUT, dpi=260, bbox_inches="tight", pad_inches=0.12)


if __name__ == "__main__":
    main()
