#!/usr/bin/env python3
"""Build the IBD overview figure for the gut-dCST paper."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_paths import FIGURES_DIR, GUT_MICROBIOME_ROOT, TABLES_DIR
from taxon_formatting import italicize_taxa_mpl

RUN_DIR = (
    GUT_MICROBIOME_ROOT
    / "outputs"
    / "dcst_analysis"
    / "runs"
    / "2026-04-11-absorb-depthscan-adaptive"
)
FOLLOWUP_DIR = RUN_DIR / "ibd_label_followup"
VALIDATION_TABLE = TABLES_DIR / "TABLE_V1_external_validation_summary.tsv"
OUT_FIG = FIGURES_DIR / "FIGURE_2_ibd_results_overview.png"

FOCAL_OUTCOMES = [
    ("IBD", "IBD", "#8c1c13"),
    ("Autoimmune", "Autoimmune disease", "#355070"),
    ("IBS", "IBS", "#6a994e"),
]

IBD_STATE_ROWS = [
    ("Bacteroides", 1, "Bacteroides"),
    ("Bacteroides__Lachnospiraceae", 2, "Bacteroides / Lachnospiraceae"),
    (
        "Bacteroides__Lachnospiraceae__Alistipes",
        3,
        "Bacteroides / Lachnospiraceae / Alistipes",
    ),
    (
        "Bacteroides__Lachnospiraceae__Alistipes__Faecalibacterium",
        4,
        "Bacteroides / Lachnospiraceae / Alistipes / Faecalibacterium",
    ),
    ("Enterobacterales", 1, "Enterobacterales"),
    ("Morganella", 1, "Morganella"),
    ("Proteus", 1, "Proteus"),
]


def load_omnibus() -> pd.DataFrame:
    omnibus = pd.read_csv(RUN_DIR / "omnibus_by_depth.tsv", sep="\t")
    return omnibus[omnibus["depth"].between(1, 4)].copy()


def load_ibd_followup() -> pd.DataFrame:
    return pd.read_csv(FOLLOWUP_DIR / "ibd_label_level_results.tsv", sep="\t")


def load_validation() -> pd.DataFrame:
    if not VALIDATION_TABLE.exists():
        raise FileNotFoundError(
            f"Missing {VALIDATION_TABLE}. Run build_validation_assets.py first."
        )
    table = pd.read_csv(VALIDATION_TABLE, sep="\t")
    for col in ("rebuilt_best_q", "frozen_best_q"):
        table[col] = pd.to_numeric(table[col].replace("NA", np.nan), errors="coerce")
    return table


def q_to_score(value: float) -> float:
    if pd.isna(value):
        return 0.0
    return -math.log10(max(float(value), 1e-12))


def pretty_label(label: str) -> str:
    return italicize_taxa_mpl(label)


def wrap_state_label(label: str, max_parts_per_line: int = 2) -> str:
    parts = [part.strip() for part in label.split(" / ")]
    chunks = [
        " / ".join(parts[idx : idx + max_parts_per_line])
        for idx in range(0, len(parts), max_parts_per_line)
    ]
    return "\n".join(pretty_label(chunk) for chunk in chunks)


def support_label(row: pd.Series, depth: int) -> str:
    cases = int(row["n_cases_in_label_all"])
    controls = int(row["n_controls_in_label_all"])
    return f"depth {depth}; IBD/control={cases}/{controls}"


def get_followup_row(df: pd.DataFrame, label: str, depth: int) -> pd.Series:
    sub = df[(df["label"] == label) & (df["depth"] == depth)]
    if sub.empty:
        raise KeyError(f"Missing IBD follow-up row for {label} at depth {depth}")
    return sub.iloc[0]


def build_figure() -> None:
    omnibus = load_omnibus()
    followup = load_ibd_followup()
    validation = load_validation()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    fig = plt.figure(figsize=(14.4, 9.3))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.0, 1.16],
        width_ratios=[0.82, 1.18],
        hspace=0.36,
        wspace=0.62,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_c = fig.add_subplot(gs[1, :])

    for key, display, color in FOCAL_OUTCOMES:
        sub = omnibus[omnibus["Condition"] == key].sort_values("depth")
        ax_a.plot(
            sub["depth"],
            sub["cramers_v"],
            marker="o",
            linewidth=2.4,
            ms=5.5,
            color=color,
            label=display,
        )
    ax_a.set_xticks([1, 2, 3, 4])
    ax_a.set_xlabel("dCST hierarchy depth")
    ax_a.set_ylabel("Cramer's V")
    ax_a.set_ylim(0.05, 0.38)
    ax_a.set_title("A. IBD is the strongest depth-resolved signal", loc="left", fontweight="bold")
    ax_a.grid(axis="y", alpha=0.22, linewidth=0.6)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.legend(frameon=False, loc="upper left")

    rows = [get_followup_row(followup, label, depth) for label, depth, _ in IBD_STATE_ROWS]
    y = np.arange(len(rows))[::-1]
    labels = [
        wrap_state_label(display) + "\n" + support_label(row, depth)
        for row, (_, depth, display) in zip(rows, IBD_STATE_ROWS)
    ]
    for idx, row in enumerate(rows):
        yy = y[idx]
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
                s=34,
                zorder=3,
                label="Frequentist adjusted OR" if idx == 0 else None,
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
            s=34,
            zorder=3,
            label="Bayesian adjusted OR" if idx == 0 else None,
        )
    ax_b.axvline(1.0, color="#6c757d", linestyle="--", linewidth=1.0)
    ax_b.set_xscale("log")
    ax_b.set_xlabel("Adjusted odds ratio for IBD")
    ax_b.set_yticks(y)
    ax_b.set_yticklabels(labels)
    ax_b.tick_params(axis="y", labelsize=8.2)
    ax_b.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax_b.spines["top"].set_visible(False)
    ax_b.spines["right"].set_visible(False)
    ax_b.set_title("B. Common dominance-lineages and rare states point to IBD dysbiosis", loc="left", fontweight="bold")
    ax_b.legend(frameon=False, loc="upper right", fontsize=8.5)

    val = validation.copy()
    y2 = np.arange(len(val))[::-1]
    rebuilt_x = [q_to_score(v) for v in val["rebuilt_best_q"]]
    frozen_x = [q_to_score(v) for v in val["frozen_best_q"]]
    labels2 = [x.replace(" ", "\n", 1) for x in val["cohort"]]
    for yy, rx, fx in zip(y2, rebuilt_x, frozen_x):
        ax_c.plot([rx, fx], [yy, yy], color="#c7c7c7", linewidth=1.4, zorder=1)
    ax_c.scatter(rebuilt_x, y2, s=58, color="#4c78a8", label="Rebuilt cohort", zorder=3)
    ax_c.scatter(frozen_x, y2, s=58, color="#f58518", label="AGP-derived label transfer", zorder=3)
    ax_c.axvline(-math.log10(0.05), color="#303030", linestyle="--", linewidth=1.0)
    ax_c.set_yticks(y2)
    ax_c.set_yticklabels(labels2, fontsize=8.5)
    ax_c.set_xlabel(r"Strongest corrected validation signal ($-\log_{10} q$)")
    ax_c.set_title(
        "C. External cohorts support IBD dominance-pattern reproducibility more strongly than universal label portability",
        loc="left",
        fontweight="bold",
    )
    ax_c.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax_c.spines["top"].set_visible(False)
    ax_c.spines["right"].set_visible(False)
    ax_c.legend(frameon=False, loc="lower right")

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


if __name__ == "__main__":
    build_figure()
