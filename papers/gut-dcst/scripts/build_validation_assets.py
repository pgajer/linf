#!/usr/bin/env python3
"""Build manuscript-facing validation table and figure assets."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

from paper_paths import DCST_VALIDATION_DIR, FIGURES_DIR, TABLES_DIR

VALIDATION_ROOT = DCST_VALIDATION_DIR

HMP2_DIR = VALIDATION_ROOT / "hmp2"
PRJEB_DIR = VALIDATION_ROOT / "prjeb84421"

TABLE_TSV = TABLES_DIR / "TABLE_V1_external_validation_summary.tsv"
TABLE_MD = TABLES_DIR / "TABLE_V1_external_validation_summary.md"
FIGURE_PNG = FIGURES_DIR / "FIGURE_V1_two_cohort_external_validation.png"


def shorten_label(label: str) -> str:
    mapping = {
        "RARE_DOMINANT": "RARE_DOMINANT",
        "Bacteroides sp.": "Bacteroides",
        "Bacteroides sp.__Faecalibacterium sp.": "Bacteroides__Faecalibacterium",
        "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__": "Faecalibacterium",
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__": "Bacteroides",
        "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Subdoligranulum;s__": "Subdoligranulum",
        "d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s____d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s__": "Faecalibacterium__Bacteroides",
        "d__Bacteria;p__Bacteroidota;c__Bacteroidia;o__Bacteroidales;f__Bacteroidaceae;g__Bacteroides;s____d__Bacteria;p__Firmicutes;c__Clostridia;o__Oscillospirales;f__Ruminococcaceae;g__Faecalibacterium;s__": "Bacteroides__Faecalibacterium",
    }
    return mapping.get(label, label)


def build_table() -> pd.DataFrame:
    rows = [
        {
            "cohort": "HMP2 / IBDMDB",
            "samples_analyzed": 166,
            "phenotype_split": "45 Healthy / 81 Crohn / 40 UC",
            "ecological_caveat": "mucosal gut sites; repeated measures",
            "strongest_depth1_signal": "no q < 0.05 depth-1 replication",
            "strongest_depth2_signal": (
                "Bacteroides__Faecalibacterium depleted in "
                "IBD_vs_Healthy, Crohn_vs_Healthy, and UC_vs_Healthy"
            ),
            "interpretation": "directional / structural replication at depth 2",
        },
        {
            "cohort": "PRJEB84421 / OFGCD-FI-2025",
            "samples_analyzed": 73,
            "phenotype_split": "20 Healthy / 24 Crohn / 29 OFG",
            "ecological_caveat": "stool-based but pediatric and phenotype-mixed",
            "strongest_depth1_signal": "RARE_DOMINANT enriched in OFG_vs_Healthy",
            "strongest_depth2_signal": "RARE_DOMINANT enriched in inflammatory_vs_healthy",
            "interpretation": "stool-based structural complement with modest cohort-specific signal",
        },
    ]
    return pd.DataFrame(rows)


def write_table(table: pd.DataFrame) -> None:
    table.to_csv(TABLE_TSV, sep="\t", index=False)
    headers = list(table.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in table.iterrows():
        vals = [str(row[h]) for h in headers]
        lines.append("| " + " | ".join(vals) + " |")
    TABLE_MD.write_text(
        "# Table V1: External validation cohort summary\n\n" + "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def add_bar_panel(ax, counts: pd.Series, title: str, color: str) -> None:
    counts = counts.sort_values(ascending=False)
    labels = [shorten_label(x) for x in counts.index]
    ax.bar(range(len(counts)), counts.values, color=color, edgecolor="black", linewidth=0.6)
    ax.set_title(title, fontsize=11)
    ax.set_ylabel("Samples")
    ax.set_xticks(range(len(counts)))
    ax.set_xticklabels(labels, rotation=35, ha="right", fontsize=8)
    ax.grid(axis="y", alpha=0.25)


def build_figure() -> None:
    hmp2_assign = pd.read_csv(HMP2_DIR / "hmp2_dcst_assignments.csv")
    prjeb_assign = pd.read_csv(PRJEB_DIR / "prjeb84421_dcst_assignments.csv")
    hmp2_d2 = pd.read_csv(HMP2_DIR / "hmp2_depth2_results.csv")
    prjeb_d1 = pd.read_csv(PRJEB_DIR / "prjeb84421_depth1_results.csv")
    prjeb_d2 = pd.read_csv(PRJEB_DIR / "prjeb84421_depth2_results.csv")

    hmp2_depth1_counts = hmp2_assign["dcst_depth1"].value_counts()
    hmp2_depth2_counts = hmp2_assign["dcst_depth2"].value_counts()
    prjeb_depth1_counts = prjeb_assign["depth1"].value_counts()

    plot_rows = [
        {
            "cohort": "HMP2",
            "label": "Bacteroides__Faecalibacterium\nIBD vs Healthy",
            "or": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "IBD_vs_Healthy"),
                "OR",
            ].iloc[0]),
            "lo": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "IBD_vs_Healthy"),
                "CI_low",
            ].iloc[0]),
            "hi": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "IBD_vs_Healthy"),
                "CI_high",
            ].iloc[0]),
            "q": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "IBD_vs_Healthy"),
                "q_value",
            ].iloc[0]),
        },
        {
            "cohort": "HMP2",
            "label": "Bacteroides__Faecalibacterium\nCrohn vs Healthy",
            "or": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "Crohn_vs_Healthy"),
                "OR",
            ].iloc[0]),
            "lo": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "Crohn_vs_Healthy"),
                "CI_low",
            ].iloc[0]),
            "hi": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "Crohn_vs_Healthy"),
                "CI_high",
            ].iloc[0]),
            "q": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "Crohn_vs_Healthy"),
                "q_value",
            ].iloc[0]),
        },
        {
            "cohort": "HMP2",
            "label": "Bacteroides__Faecalibacterium\nUC vs Healthy",
            "or": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "UC_vs_Healthy"),
                "OR",
            ].iloc[0]),
            "lo": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "UC_vs_Healthy"),
                "CI_low",
            ].iloc[0]),
            "hi": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "UC_vs_Healthy"),
                "CI_high",
            ].iloc[0]),
            "q": float(hmp2_d2.loc[
                (hmp2_d2["DCST"] == "Bacteroides sp.__Faecalibacterium sp.")
                & (hmp2_d2["Condition"] == "UC_vs_Healthy"),
                "q_value",
            ].iloc[0]),
        },
        {
            "cohort": "PRJEB84421",
            "label": "RARE_DOMINANT\nInflammatory vs Healthy",
            "or": float(prjeb_d2.loc[
                (prjeb_d2["DCST"] == "RARE_DOMINANT")
                & (prjeb_d2["Condition"] == "Inflammatory_vs_Healthy"),
                "OR",
            ].iloc[0]),
            "lo": float(prjeb_d2.loc[
                (prjeb_d2["DCST"] == "RARE_DOMINANT")
                & (prjeb_d2["Condition"] == "Inflammatory_vs_Healthy"),
                "CI_low",
            ].iloc[0]),
            "hi": float(prjeb_d2.loc[
                (prjeb_d2["DCST"] == "RARE_DOMINANT")
                & (prjeb_d2["Condition"] == "Inflammatory_vs_Healthy"),
                "CI_high",
            ].iloc[0]),
            "q": float(prjeb_d2.loc[
                (prjeb_d2["DCST"] == "RARE_DOMINANT")
                & (prjeb_d2["Condition"] == "Inflammatory_vs_Healthy"),
                "q_value",
            ].iloc[0]),
        },
        {
            "cohort": "PRJEB84421",
            "label": "RARE_DOMINANT\nOFG vs Healthy",
            "or": float(prjeb_d1.loc[
                (prjeb_d1["DCST"] == "RARE_DOMINANT")
                & (prjeb_d1["Condition"] == "OFG_vs_Healthy"),
                "OR",
            ].iloc[0]),
            "lo": float(prjeb_d1.loc[
                (prjeb_d1["DCST"] == "RARE_DOMINANT")
                & (prjeb_d1["Condition"] == "OFG_vs_Healthy"),
                "CI_low",
            ].iloc[0]),
            "hi": float(prjeb_d1.loc[
                (prjeb_d1["DCST"] == "RARE_DOMINANT")
                & (prjeb_d1["Condition"] == "OFG_vs_Healthy"),
                "CI_high",
            ].iloc[0]),
            "q": float(prjeb_d1.loc[
                (prjeb_d1["DCST"] == "RARE_DOMINANT")
                & (prjeb_d1["Condition"] == "OFG_vs_Healthy"),
                "q_value",
            ].iloc[0]),
        },
    ]
    plot_df = pd.DataFrame(plot_rows)

    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.15], hspace=0.45, wspace=0.28)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax3 = fig.add_subplot(gs[1, 0])
    ax4 = fig.add_subplot(gs[1, 1])

    add_bar_panel(ax1, hmp2_depth1_counts, "A. HMP2 Depth-1 DCST Sizes", "#4C78A8")
    add_bar_panel(ax2, hmp2_depth2_counts.head(6), "B. HMP2 Depth-2 DCST Sizes", "#72B7B2")
    add_bar_panel(ax3, prjeb_depth1_counts.head(6), "C. PRJEB84421 Depth-1 DCST Sizes", "#F28E2B")

    y = list(range(len(plot_df)))[::-1]
    colors = ["#4C78A8" if c == "HMP2" else "#F28E2B" for c in plot_df["cohort"]]
    ax4.axvline(1.0, color="black", linestyle="--", linewidth=1)
    for yi, (_, row) in zip(y, plot_df.iterrows()):
        ax4.plot([row["lo"], row["hi"]], [yi, yi], color="gray", linewidth=1.5)
        ax4.scatter(row["or"], yi, s=55, color=colors[len(y) - 1 - yi], zorder=3)
        ax4.text(
            row["hi"] * 1.06,
            yi,
            f"q={row['q']:.3f}",
            va="center",
            fontsize=8,
        )
    ax4.set_yticks(y)
    ax4.set_yticklabels(plot_df["label"], fontsize=8)
    ax4.set_xscale("log")
    ax4.set_xlabel("Odds ratio (log scale)")
    ax4.set_title("D. Directional Replication Highlights", fontsize=11)
    ax4.grid(axis="x", alpha=0.25)

    fig.suptitle(
        "Figure V1. Two-cohort external validation overview",
        fontsize=14,
        y=0.98,
    )
    fig.text(
        0.5,
        0.015,
        "HMP2 provides clinically grounded IBD validation with strongest signals at depth 2; "
        "PRJEB84421 provides a stool-based inflammatory complement with more modest cohort-specific signals.",
        ha="center",
        fontsize=9,
    )
    fig.savefig(FIGURE_PNG, dpi=200, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    table = build_table()
    write_table(table)
    build_figure()
    print(f"Wrote {TABLE_TSV}")
    print(f"Wrote {TABLE_MD}")
    print(f"Wrote {FIGURE_PNG}")


if __name__ == "__main__":
    main()
