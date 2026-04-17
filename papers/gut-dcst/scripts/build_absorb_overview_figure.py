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
)
OUT_FIG = FIGURES_DIR / "FIGURE_2_absorb_overview.png"

FOCAL_OUTCOMES = [
    ("IBD", "IBD", "#8c1c13"),
    ("Autoimmune", "Autoimmune", "#355070"),
    ("IBS", "IBS", "#6a994e"),
]


def pretty_label(label: str) -> str:
    return italicize_taxa_mpl(label.replace("__", " / "))


def load_top_depth1() -> tuple[list[str], list[int]]:
    df = pd.read_csv(RUN_DIR / "largest_labels_by_depth.tsv", sep="\t")
    sub = df[(df["depth"] == 1) & (df["rank"] <= 10)].copy()
    labels = [pretty_label(x) for x in sub["label"]]
    counts = sub["n"].astype(int).tolist()
    return labels, counts


def load_depth_summary() -> pd.DataFrame:
    return pd.read_csv(RUN_DIR / "depth_summary.tsv", sep="\t")


def load_omnibus() -> pd.DataFrame:
    omni = pd.read_csv(RUN_DIR / "omnibus_by_depth.tsv", sep="\t")
    return omni[omni["depth"].between(1, 4)].copy()


def build_figure() -> None:
    top_labels, top_counts = load_top_depth1()
    depth_summary = load_depth_summary()
    omnibus = load_omnibus()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    fig = plt.figure(figsize=(13.4, 8.9))
    gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.05], hspace=0.36, wspace=0.28)

    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    colors = [
        "#355070",
        "#588157",
        "#bc6c25",
        "#8c1c13",
        "#6d597a",
        "#457b9d",
        "#adb5bd",
        "#adb5bd",
        "#adb5bd",
        "#adb5bd",
    ]
    y = np.arange(len(top_labels))
    bars = ax_a.barh(y, top_counts, color=colors[: len(top_labels)], edgecolor="#2f2f2f", linewidth=0.7)
    ax_a.set_yticks(y)
    ax_a.set_yticklabels(top_labels)
    ax_a.invert_yaxis()
    ax_a.set_xlabel("Samples")
    ax_a.set_title("A. Dominant absorb states at depth 1", loc="left", fontweight="bold")
    ax_a.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    for bar, count in zip(bars[:6], top_counts[:6]):
        ax_a.text(
            bar.get_width() + 120,
            bar.get_y() + bar.get_height() / 2,
            f"{count:,}",
            va="center",
            ha="left",
            fontsize=8.8,
            color="#303030",
        )
    ax_a.text(
        0.995,
        1.03,
        "30,290 samples; 40 named absorb states",
        transform=ax_a.transAxes,
        ha="right",
        va="bottom",
        fontsize=9.5,
        color="#4f4f4f",
    )

    depths = depth_summary["depth"].astype(int).to_numpy()
    n_labels = depth_summary["n_labels"].to_numpy()
    share_changed = depth_summary["share_changed_from_rare_view"].to_numpy() * 100.0
    ax_b.plot(depths, n_labels, marker="o", color="#355070", linewidth=2.4, ms=6)
    ax_b.set_xlabel("Hierarchy depth")
    ax_b.set_ylabel("Named absorb states", color="#355070")
    ax_b.tick_params(axis="y", labelcolor="#355070")
    ax_b.set_xticks(depths)
    ax_b.set_title("B. Hierarchy growth and divergence from the pure view", loc="left", fontweight="bold")
    ax_b.grid(axis="y", alpha=0.22, linewidth=0.6)
    ax_b.spines["top"].set_visible(False)
    ax_b2 = ax_b.twinx()
    ax_b2.plot(depths, share_changed, marker="s", color="#bc6c25", linewidth=2.1, ms=5)
    ax_b2.set_ylabel("Assignments changed vs pure hierarchy (%)", color="#bc6c25")
    ax_b2.tick_params(axis="y", labelcolor="#bc6c25")
    ax_b2.spines["top"].set_visible(False)
    ax_b.axvspan(0.85, 4.15, color="#d9ead3", alpha=0.18)
    ax_b.text(
        0.03,
        0.08,
        "Primary manuscript depths",
        transform=ax_b.transAxes,
        fontsize=9,
        color="#4f4f4f",
    )
    for x, yv in zip(depths[:4], n_labels[:4]):
        ax_b.text(x, yv + 18, f"{yv}", ha="center", va="bottom", fontsize=8.5, color="#355070")

    focal_depths = [1, 2, 3, 4]
    for key, display, color in FOCAL_OUTCOMES:
        sub = omnibus[omnibus["Condition"] == key].sort_values("depth")
        ax_c.plot(
            sub["depth"],
            sub["cramers_v"],
            marker="o",
            linewidth=2.2,
            ms=5,
            color=color,
            label=display,
        )
    ax_c.set_xticks(focal_depths)
    ax_c.set_xlabel("Hierarchy depth")
    ax_c.set_ylabel("Cramer's V")
    ax_c.set_ylim(0.05, 0.38)
    ax_c.set_title("C. Absorb-state association strengthens in the main inflammatory branches", loc="left", fontweight="bold")
    ax_c.grid(axis="y", alpha=0.22, linewidth=0.6)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.legend(ncol=3, frameon=False, loc="upper left")

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220, bbox_inches="tight", pad_inches=0.04)
    plt.close(fig)


if __name__ == "__main__":
    build_figure()
