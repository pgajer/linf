from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_paths import FIGURES_DIR, GUT_MICROBIOME_ROOT
from taxon_formatting import italicize_taxa_mpl

RUN_DIR = (
    GUT_MICROBIOME_ROOT
    / "outputs"
    / "dcst_analysis"
    / "runs"
    / "2026-04-11-absorb-depthscan-adaptive"
    / "ibd_label_followup"
)
OUT_FIG = FIGURES_DIR / "FIGURE_3_ibd_absorb_dominance_lineages.png"


def load_results() -> pd.DataFrame:
    return pd.read_csv(RUN_DIR / "ibd_label_level_results.tsv", sep="\t")


def get_row(df: pd.DataFrame, label: str, depth: int) -> pd.Series:
    sub = df[(df["label"] == label) & (df["depth"] == depth)]
    if sub.empty:
        raise KeyError(f"Missing IBD row for {label} at depth {depth}")
    return sub.iloc[0]


def pretty_label(label: str) -> str:
    return italicize_taxa_mpl(label.replace("__", " / "))


def wrap_lineage_label(label: str, max_parts_per_line: int = 2) -> str:
    parts = label.split("__")
    chunks = []
    for idx in range(0, len(parts), max_parts_per_line):
        chunk = " / ".join(parts[idx : idx + max_parts_per_line])
        chunks.append(pretty_label(chunk))
    return "\n".join(chunks)


def support_label(row: pd.Series) -> str:
    cases = int(row["n_cases_in_label_all"])
    controls = int(row["n_controls_in_label_all"])
    total = int(row["n_label_all"])
    return f"n={total}; IBD/control={cases}/{controls}"


def build_figure() -> None:
    df = load_results()

    lineage_rows = [
        get_row(df, "Bacteroides", 1),
        get_row(df, "Bacteroides__Lachnospiraceae", 2),
        get_row(df, "Bacteroides__Lachnospiraceae__Alistipes", 3),
        get_row(df, "Bacteroides__Lachnospiraceae__Alistipes__Faecalibacterium", 4),
    ]

    forest_rows = [
        get_row(df, "Bacteroides", 1),
        get_row(df, "Enterobacterales", 1),
        get_row(df, "Morganella", 1),
        get_row(df, "Proteus", 1),
        get_row(df, "Bacteroides__Lachnospiraceae", 2),
        get_row(df, "Bacteroides__Lachnospiraceae__Alistipes", 3),
        get_row(df, "Bacteroides__Lachnospiraceae__Alistipes__Faecalibacterium", 4),
    ]

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    fig = plt.figure(figsize=(13.6, 8.8))
    gs = fig.add_gridspec(1, 2, width_ratios=[0.88, 1.52], wspace=0.55)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    ax_a.set_title("A. IBD signal resolves into\na deeper absorb dominance-lineage", loc="left", fontweight="bold")
    ax_a.set_axis_off()

    box_colors = ["#d9ead3", "#cfe2f3", "#fce5cd", "#ead1dc"]
    y_positions = [0.88, 0.65, 0.42, 0.19]
    x_center = 0.5
    box_w = 0.90
    box_h = 0.14

    for idx, (row, y, face) in enumerate(zip(lineage_rows, y_positions, box_colors)):
        label = wrap_lineage_label(row["label"])
        ax_a.add_patch(
            plt.Rectangle(
                (x_center - box_w / 2, y - box_h / 2),
                box_w,
                box_h,
                facecolor=face,
                edgecolor="#4f4f4f",
                linewidth=1.0,
            )
        )
        ax_a.text(
            x_center,
            y + 0.027,
            label,
            ha="center",
            va="center",
            fontsize=10.0,
            fontweight="bold" if idx == 0 else None,
        )
        ax_a.text(
            x_center,
            y - 0.042,
            f"Depth {int(row['depth'])} | n = {int(row['n_label_all'])} | adjusted OR = {row['freq_adj_or']:.2f}",
            ha="center",
            va="center",
            fontsize=8.9,
            color="#333333",
        )
        if idx < len(lineage_rows) - 1:
            ax_a.annotate(
                "",
                xy=(x_center, y_positions[idx + 1] + box_h / 2),
                xytext=(x_center, y - box_h / 2),
                arrowprops=dict(arrowstyle="-|>", color="#5f5f5f", linewidth=1.5),
            )
    ax_b.set_title("B. Counts and adjusted estimates\nfor representative IBD states", loc="left", fontweight="bold")
    y = np.arange(len(forest_rows))[::-1]
    labels = [wrap_lineage_label(r["label"], max_parts_per_line=1) + "\n" + support_label(r) for r in forest_rows]

    for i, row in enumerate(forest_rows):
        yy = y[i]
        if np.isfinite(row["freq_adj_or"]):
            ax_b.plot(
                [row["freq_adj_ci_lo"], row["freq_adj_ci_hi"]],
                [yy + 0.12, yy + 0.12],
                color="#355070",
                linewidth=2.0,
            )
            ax_b.scatter(
                row["freq_adj_or"],
                yy + 0.12,
                color="#355070",
                s=38,
                zorder=3,
                label="Frequentist adjusted OR" if i == 0 else None,
            )
        ax_b.plot(
            [row["bayes_adj_ci_lo"], row["bayes_adj_ci_hi"]],
            [yy - 0.12, yy - 0.12],
            color="#bc6c25",
            linewidth=2.0,
        )
        ax_b.scatter(
            row["bayes_adj_or"],
            yy - 0.12,
            color="#bc6c25",
            s=38,
            zorder=3,
            label="Bayesian adjusted OR" if i == 0 else None,
        )

    ax_b.axvline(1.0, color="#6c757d", linestyle="--", linewidth=1.1)
    ax_b.set_xscale("log")
    ax_b.set_xlabel("Adjusted odds ratio")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(labels)
    ax_b.tick_params(axis="y", labelsize=8.6)
    ax_b.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.legend(frameon=False, loc="upper right")
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    build_figure()
