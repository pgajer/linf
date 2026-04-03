from __future__ import annotations

import csv
import gzip
from collections import Counter
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_paths import DCST_ANALYSIS_DIR, FIGURES_DIR, GUT_MICROBIOME_ROOT

META_FILE = GUT_MICROBIOME_ROOT / "data" / "prime_gut_project_sample_metadata_2026-03-24.csv.gz"
ABUND_FILE = GUT_MICROBIOME_ROOT / "outputs" / "prime_species" / "prime_gut_projects_silva_species_absolute_2026-03-24.csv.gz"
ASSIGN_FILE = DCST_ANALYSIS_DIR / "full_cohort_dcst_assignments.csv"
OUT_FIG = FIGURES_DIR / "FIGURE_2_full_cohort_landscape.png"

AGP_PROJECT = "PRJEB11419"
MIN_LIB = 1000
PREV_PROP = 0.05
MIN_COUNT = 2


def load_agp_runs() -> set[str]:
    meta = pd.read_csv(META_FILE, compression="gzip", low_memory=False, usecols=["Run", "BioProject"])
    return set(meta.loc[meta["BioProject"] == AGP_PROJECT, "Run"].astype(str))


def load_agp_counts(agp_runs: set[str]) -> np.ndarray:
    rows: list[list[int]] = []
    with gzip.open(ABUND_FILE, "rt") as f:
        reader = csv.reader(f)
        header = next(reader)
        n_taxa = len(header) - 1
        for row in reader:
            if row[0] in agp_runs:
                rows.append([int(x) for x in row[1:]])
    counts = np.asarray(rows, dtype=np.int64)
    if counts.ndim != 2 or counts.shape[1] != n_taxa:
        raise RuntimeError("Unexpected AGP count matrix shape")
    return counts


def filter_counts(counts: np.ndarray) -> np.ndarray:
    lib_sizes = counts.sum(axis=1)
    counts = counts[lib_sizes >= MIN_LIB]
    present = counts >= MIN_COUNT
    prev_threshold = int(np.ceil(PREV_PROP * counts.shape[0]))
    keep_feat = present.sum(axis=0) >= prev_threshold
    counts = counts[:, keep_feat]
    counts = counts[counts.sum(axis=1) > 0]
    return counts


def short_label(label: str) -> str:
    return label.replace("d__", "")


def load_depth1_sizes() -> tuple[list[str], list[int]]:
    df = pd.read_csv(ASSIGN_FILE, usecols=["dcst_depth1_short"])
    size_table = Counter(df["dcst_depth1_short"].astype(str))
    ordered = sorted(size_table.items(), key=lambda kv: (-kv[1], kv[0]))
    labels = [short_label(k) for k, _ in ordered]
    counts = [v for _, v in ordered]
    return labels, counts


def build_figure(labels: list[str], size_counts: list[int], dom_values: np.ndarray) -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    fig = plt.figure(figsize=(13, 8.6))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.34)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    top_n = 12
    tail_count = int(sum(size_counts[top_n:]))
    plot_labels = labels[:top_n] + ["Other named DCSTs"]
    plot_counts = size_counts[:top_n] + [tail_count]
    colors = ["#335c67", "#335c67", "#7cb518", "#ff9f1c", "#c8553d"] + ["#7a8c99"] * (len(plot_labels) - 5)
    if len(colors) < len(plot_labels):
        colors += ["#7a8c99"] * (len(plot_labels) - len(colors))
    colors[-1] = "#adb5bd"

    y = np.arange(len(plot_labels))
    bars = ax1.barh(y, plot_counts, color=colors, edgecolor="#2f2f2f", linewidth=0.7)
    ax1.set_title("A. Depth-1 DCST landscape in the AGP full cohort", loc="left", fontweight="bold")
    ax1.set_xlabel("Samples")
    ax1.set_yticks(y)
    ax1.set_yticklabels(plot_labels)
    ax1.invert_yaxis()
    ax1.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.text(
        0.995,
        1.03,
        "30,290 samples; 40 named DCSTs + RARE_DOMINANT",
        transform=ax1.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        color="#4d4d4d",
    )
    for idx in range(min(6, len(plot_labels))):
        ax1.text(
            bars[idx].get_width() + 110,
            bars[idx].get_y() + bars[idx].get_height() / 2,
            f"{plot_counts[idx]:,}",
            ha="left",
            va="center",
            fontsize=8.8,
            color="#2f2f2f",
        )

    ax2.hist(dom_values, bins=50, color="#5abf90", edgecolor="white", linewidth=0.8)
    ax2.axvline(0.5, color="#c8553d", linestyle="--", linewidth=1.5)
    ax2.text(0.505, ax2.get_ylim()[1] * 0.93, "50% dominance", color="#c8553d", fontsize=9, ha="left")
    q25, q50, q75 = np.quantile(dom_values, [0.25, 0.5, 0.75])
    for val, label, col in [(q25, "Q1", "#4d908e"), (q50, "Median", "#1d3557"), (q75, "Q3", "#4d908e")]:
        ax2.axvline(val, color=col, linestyle=":", linewidth=1.3)
        ax2.text(val + 0.01, ax2.get_ylim()[1] * 0.82, label, color=col, fontsize=8.5)
    ax2.set_title("B. Empirical dominance strength in the AGP cohort", loc="left", fontweight="bold")
    ax2.set_xlabel("Relative abundance of dominant species")
    ax2.set_ylabel("Samples")
    ax2.grid(axis="y", alpha=0.25, linewidth=0.6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.text(
        0.995,
        1.03,
        f"median = {q50:.2f}; IQR = [{q25:.2f}, {q75:.2f}]",
        transform=ax2.transAxes,
        ha="right",
        va="bottom",
        fontsize=10,
        color="#4d4d4d",
    )

    fig.suptitle(
        "Figure 2. The full-cohort AGP gut landscape combines a few large dominant states with a broad moderate-dominance background",
        fontsize=14.5,
        fontweight="bold",
        y=0.98,
    )
    fig.text(
        0.5,
        0.015,
        "Panel A shows the leading depth-1 DCSTs after truncation at n0 = 50. Panel B shows the empirical dominance strength, defined as the relative abundance of the most abundant taxon in each sample.",
        ha="center",
        va="bottom",
        fontsize=10.2,
        color="#3b3b3b",
    )
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    agp_runs = load_agp_runs()
    counts = load_agp_counts(agp_runs)
    counts = filter_counts(counts)
    if counts.shape != (30290, 274):
        print(f"Warning: filtered matrix shape is {counts.shape}, expected about (30290, 274)")
    dom_values = counts.max(axis=1) / counts.sum(axis=1)
    labels, size_counts = load_depth1_sizes()
    build_figure(labels, size_counts, dom_values)


if __name__ == "__main__":
    main()
