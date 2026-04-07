#!/usr/bin/env python3
"""Build a manuscript-ready four-cohort external validation figure."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_paths import DCST_VALIDATION_DIR, FIGURES_DIR
from taxon_formatting import italicize_taxa_mpl

REPORT_DIR = DCST_VALIDATION_DIR / "completed_validation_report_assets_2026-04-02"
OUT = FIGURES_DIR / "FIGURE_V1_four_cohort_external_validation.png"

COHORT_ORDER = ["hmp2", "prjeb84421", "halfvarson_2017", "gevers_2014"]
DISPLAY = {
    "hmp2": "HMP2 / IBDMDB",
    "prjeb84421": "PRJEB84421 / OFGCD-FI-2025",
    "halfvarson_2017": "Halfvarson 2017",
    "gevers_2014": "Gevers 2014",
}
PRIMARY = {
    "hmp2": "IBD_vs_Healthy",
    "prjeb84421": "Inflammatory_vs_Healthy",
    "halfvarson_2017": "IBD_vs_Healthy",
    "gevers_2014": "Crohn_vs_Healthy",
}
ROLE = {
    "hmp2": "Clinically grounded\nIBD validation",
    "prjeb84421": "Stool inflammatory\ncomplement",
    "halfvarson_2017": "Strong external\npositive validation",
    "gevers_2014": "Completed external\nweak/null control",
}


def compact_state(label: str) -> str:
    if label == "RARE_DOMINANT":
        return "Rare"
    label = label.replace("Bacteroides sp.__Faecalibacterium sp.", "Bacteroides/Faecalibacterium")
    label = label.replace(
        "Segatella sp.__Faecalibacterium prausnitzii",
        "Segatella/F. prausnitzii",
    )
    label = label.replace(
        "Faecalibacterium prausnitzii__Bacteroides vulgatus",
        "F. prausnitzii/B. vulgatus",
    )
    label = label.replace("Prevotella_9 sp.__Prevotella_7 sp.", "Prevotella9/Prevotella7")
    label = label.replace("Bacteroides vulgatus", "B. vulgatus")
    label = label.replace("Segatella sp.", "Segatella")
    label = label.replace("Prevotella_9 sp.", "Prevotella9")
    label = label.replace("Subdoligranulum sp.", "Subdoligranulum")
    label = label.replace("Bacteroides sp.__Alistipes sp.", "Bacteroides/Alistipes")
    if len(label) > 28:
        return label[:25] + "..."
    return label


def main() -> None:
    cohort_df = pd.read_csv(REPORT_DIR / "cohort_overview.tsv", sep="\t")
    comp_df = pd.read_csv(REPORT_DIR / "comparison_signals.tsv", sep="\t")

    cohort_df = cohort_df.set_index("cohort_id").loc[COHORT_ORDER].reset_index()
    rows = []
    for cohort_id in COHORT_ORDER:
        row = comp_df[
            (comp_df["cohort_id"] == cohort_id)
            & (comp_df["comparison"] == PRIMARY[cohort_id])
        ].iloc[0]
        rows.append(
            {
                "cohort_id": cohort_id,
                "display_name": DISPLAY[cohort_id],
                "samples_after_qc": int(cohort_df.loc[cohort_df["cohort_id"] == cohort_id, "samples_after_qc"].iloc[0]),
                "depth1_q": float(row["depth1_best_q"]),
                "depth2_q": float(row["depth2_best_q"]),
                "depth1_state": compact_state(str(row["depth1_best_dcst"])),
                "depth2_state": compact_state(str(row["depth2_best_dcst"])),
                "depth1_hits": int(row["depth1_sig_hits"]),
                "depth2_hits": int(row["depth2_sig_hits"]),
                "role": ROLE[cohort_id],
            }
        )
    plot_df = pd.DataFrame(rows)
    plot_df["depth1_score"] = -np.log10(plot_df["depth1_q"].clip(lower=1e-12))
    plot_df["depth2_score"] = -np.log10(plot_df["depth2_q"].clip(lower=1e-12))

    fig = plt.figure(figsize=(11.6, 7.2))
    gs = fig.add_gridspec(
        2,
        2,
        height_ratios=[1.18, 0.92],
        width_ratios=[1.02, 1.18],
        hspace=0.32,
        wspace=0.32,
    )

    ax1 = fig.add_subplot(gs[0, 0])
    y = np.arange(len(plot_df))
    ax1.barh(y, plot_df["samples_after_qc"], color="#4c78a8")
    ax1.set_yticks(y, labels=plot_df["display_name"])
    ax1.invert_yaxis()
    ax1.set_xlabel("Samples after QC")
    ax1.set_title("A. Validation cohort sizes", loc="left", fontweight="bold")
    for yi, val in zip(y, plot_df["samples_after_qc"]):
        ax1.text(val + 3, yi, str(val), va="center", fontsize=8)
    ax1.set_xlim(0, plot_df["samples_after_qc"].max() + 28)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)

    ax2 = fig.add_subplot(gs[0, 1])
    ax2.scatter(plot_df["depth1_score"], y, color="#4c78a8", s=55, label="Depth 1")
    ax2.scatter(plot_df["depth2_score"], y, color="#e45756", marker="s", s=55, label="Depth 2")
    ax2.axvline(-math.log10(0.05), color="#666666", linestyle="--", linewidth=1, label="q = 0.05")
    ax2.set_yticks(y, labels=plot_df["display_name"])
    ax2.invert_yaxis()
    ax2.set_xlabel(r"$-\log_{10}(q)$ for primary contrast")
    ax2.set_title("B. Primary-comparison support", loc="left", fontweight="bold")
    ax2.legend(loc="lower right", fontsize=8, frameon=True, borderpad=0.35)
    for idx, row in plot_df.iterrows():
        ax2.text(
            row["depth1_score"] + 0.07,
            idx - 0.14,
            italicize_taxa_mpl(row["depth1_state"]),
            fontsize=7,
            color="#4c78a8",
        )
        ax2.text(
            row["depth2_score"] + 0.07,
            idx + 0.18,
            italicize_taxa_mpl(row["depth2_state"]),
            fontsize=7,
            color="#e45756",
        )
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.set_xlim(-0.35, max(plot_df["depth1_score"].max(), plot_df["depth2_score"].max()) + 0.5)

    ax3 = fig.add_subplot(gs[1, :])
    ax3.axis("off")
    ax3.set_title("C. Cohort roles and signal counts", loc="left", fontweight="bold", pad=6)
    box_specs = [
        (0.03, 0.53),
        (0.53, 0.53),
        (0.03, 0.09),
        (0.53, 0.09),
    ]
    for (x, y0), (_, row) in zip(box_specs, plot_df.iterrows()):
        text = (
            f"{row['display_name']}\n"
            f"{row['role']}\n"
            f"D1 hits: {row['depth1_hits']}   D2 hits: {row['depth2_hits']}"
        )
        ax3.text(
            x,
            y0,
            text,
            ha="left",
            va="bottom",
            fontsize=10,
            linespacing=1.18,
            bbox=dict(boxstyle="round,pad=0.55", facecolor="#f7f7f7", edgecolor="#999999"),
            transform=ax3.transAxes,
        )

    fig.subplots_adjust(left=0.08, right=0.985, top=0.96, bottom=0.08)
    fig.savefig(OUT, dpi=220, bbox_inches="tight", pad_inches=0.03)


if __name__ == "__main__":
    main()
