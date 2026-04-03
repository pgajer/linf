from __future__ import annotations

from pathlib import Path
import re

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_paths import DCST_ANALYSIS_DIR, FIGURES_DIR, PHASE1_ARCHIVE_DIR

BASE = DCST_ANALYSIS_DIR
SILVA = BASE / "full_cohort_dcst_assignments.csv"
GG2 = BASE / "gg2_dcst_assignments.csv"
CROSSTAB = BASE / "silva_vs_gg2_crosstab.csv"
SUMMARY = PHASE1_ARCHIVE_DIR / "phase1_summary_report.txt"
OUT = FIGURES_DIR / "FIGURE_6_taxonomy_concordance.png"


def pretty_label(x: str) -> str:
    mapping = {
        "RARE_DOMINANT": "Rare dominant",
        "d__Bacteria;__;__;__;__;__;__": "Unresolved",
        "Unassigned;__;__;__;__;__;__": "Unassigned",
        "d__Bacteroides_H_857956": "Bacteroides_H",
        "d__Phocaeicola_A": "Phocaeicola_A",
        "d__Enterobacteriaceae_A_725029": "Enterobacteriaceae_A",
        "d__Enterobacterales_737866": "Enterobacterales",
        "d__Agathobacter_164117": "Agathobacter",
        "d__Pseudomonas_E_647464": "Pseudomonas_E",
        "d__Prevotella_9": "Prevotella_9",
    }
    if x in mapping:
        return mapping[x]
    if x.startswith("d__"):
        return x.replace("d__", "")
    return x


def load_agreement() -> tuple[int, int, int]:
    text = SUMMARY.read_text()
    m = re.search(r"Exact genus match:\s+(\d+)/(\d+)\s+\(([\d.]+)%\)", text)
    if not m:
        raise RuntimeError("Could not parse harmonized agreement from phase1_summary_report.txt")
    harmonized = int(m.group(1))
    total = int(m.group(2))

    silva = pd.read_csv(SILVA, usecols=["Run", "dcst_depth1_short"])
    gg2 = pd.read_csv(GG2, usecols=["Run", "gg2_dcst_short"])
    merged = silva.merge(gg2, on="Run", how="inner")
    exact_short = int((merged["dcst_depth1_short"] == merged["gg2_dcst_short"]).sum())
    return harmonized, exact_short, total


def build_mapping_table() -> pd.DataFrame:
    ct = pd.read_csv(CROSSTAB)
    selected_rows = [
        "d__Akkermansia",
        "d__Faecalibacterium",
        "d__Prevotella_9",
        "d__Streptococcus",
        "d__Staphylococcus",
        "d__Enterobacteriaceae",
        "d__Bacteroides",
        "d__Escherichia-Shigella",
    ]
    categories = [
        "d__Akkermansia",
        "d__Faecalibacterium",
        "d__Prevotella",
        "d__Streptococcus",
        "d__Staphylococcus",
        "d__Enterobacteriaceae_A_725029",
        "d__Enterobacterales_737866",
        "d__Phocaeicola_A",
        "d__Bacteroides_H_857956",
        "d__Bacteria;__;__;__;__;__;__",
        "RARE_DOMINANT",
    ]

    records = []
    for row_name in selected_rows:
        row = ct.loc[ct["dcst_depth1_short"] == row_name].iloc[0]
        vals = row.drop(labels=["dcst_depth1_short"]).astype(float)
        total = vals.sum()
        out = {"silva": row_name}
        captured = 0.0
        for cat in categories:
            val = float(vals.get(cat, 0.0))
            out[cat] = val / total if total else 0.0
            captured += val
        out["Other"] = max(0.0, 1.0 - captured / total) if total else 0.0
        records.append(out)
    return pd.DataFrame.from_records(records)


def build_figure(mapping: pd.DataFrame, harmonized: int, exact_short: int, total: int) -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    category_order = [
        "d__Akkermansia",
        "d__Faecalibacterium",
        "d__Prevotella",
        "d__Streptococcus",
        "d__Staphylococcus",
        "d__Enterobacteriaceae_A_725029",
        "d__Enterobacterales_737866",
        "d__Phocaeicola_A",
        "d__Bacteroides_H_857956",
        "d__Bacteria;__;__;__;__;__;__",
        "RARE_DOMINANT",
        "Other",
    ]
    colors = {
        "d__Akkermansia": "#4daf4a",
        "d__Faecalibacterium": "#66c2a5",
        "d__Prevotella": "#8da0cb",
        "d__Streptococcus": "#fc8d62",
        "d__Staphylococcus": "#e78ac3",
        "d__Enterobacteriaceae_A_725029": "#ffd92f",
        "d__Enterobacterales_737866": "#e5c494",
        "d__Phocaeicola_A": "#1b9e77",
        "d__Bacteroides_H_857956": "#7570b3",
        "d__Bacteria;__;__;__;__;__;__": "#4c78a8",
        "RARE_DOMINANT": "#9e9e9e",
        "Other": "#d9d9d9",
    }

    fig = plt.figure(figsize=(11.2, 9.2))
    gs = fig.add_gridspec(2, 1, height_ratios=[0.88, 1.75], hspace=0.42)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    summary_df = pd.DataFrame(
        {
            "metric": ["Harmonized genus agreement", "Exact short-label match"],
            "count": [harmonized, exact_short],
        }
    )
    summary_df["fraction"] = summary_df["count"] / total
    ax1.set_title("A. Cross-taxonomy agreement depends on whether relabeling is harmonized", loc="left", fontweight="bold")
    ax1.barh(
        summary_df["metric"],
        summary_df["fraction"] * 100,
        color=["#4daf4a", "#bdbdbd"],
        edgecolor="#404040",
        height=0.46,
    )
    ax1.set_xlim(0, 100)
    ax1.set_xlabel("Samples (%)")
    ax1.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    for y, (_, row) in enumerate(summary_df.iterrows()):
        ax1.text(
            row["fraction"] * 100 + 1.2,
            y,
            f"{row['count']:,} / {total:,} ({row['fraction']*100:.1f}%)",
            va="center",
            ha="left",
            fontsize=11,
        )

    plot_df = mapping.copy()
    plot_df["silva_pretty"] = plot_df["silva"].map(pretty_label)
    y = np.arange(len(plot_df))
    left = np.zeros(len(plot_df))
    ax2.set_title("B. Most dominant states map cleanly, while a few large states split under GG2", loc="left", fontweight="bold")

    for cat in category_order:
        vals = plot_df[cat].to_numpy() * 100
        ax2.barh(
            y,
            vals,
            left=left,
            color=colors[cat],
            edgecolor="white",
            height=0.72,
            label=pretty_label(cat),
        )
        for i, v in enumerate(vals):
            if v >= 12:
                ax2.text(left[i] + v / 2, i, f"{pretty_label(cat)}\n{v:.0f}%", ha="center", va="center", fontsize=8.8, color="#1f1f1f")
        left += vals

    ax2.set_yticks(y)
    ax2.set_yticklabels(plot_df["silva_pretty"])
    ax2.set_xlim(0, 100)
    ax2.set_xlabel("Within-SILVA row share mapped to GG2")
    ax2.invert_yaxis()
    ax2.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.axhline(5.5, color="#b0b0b0", linewidth=1.0, linestyle="--")
    ax2.text(101.5, 2.5, "Near one-to-one\nmappings", va="center", ha="left", fontsize=10.3, color="#4a4a4a")
    ax2.text(101.5, 6.5, "Meaningful\nrelabelings", va="center", ha="left", fontsize=10.3, color="#4a4a4a")

    handles, labels = ax2.get_legend_handles_labels()
    keep = [
        "Akkermansia",
        "Faecalibacterium",
        "Prevotella",
        "Streptococcus",
        "Staphylococcus",
        "Enterobacteriaceae_A",
        "Enterobacterales",
        "Phocaeicola_A",
        "Bacteroides_H",
        "Unresolved",
        "Rare dominant",
        "Other",
    ]
    filtered = [(h, l) for h, l in zip(handles, labels) if l in keep]
    ax2.legend(
        [h for h, _ in filtered],
        [l for _, l in filtered],
        ncol=4,
        loc="upper center",
        bbox_to_anchor=(0.5, -0.18),
        frameon=False,
        fontsize=9.5,
        columnspacing=1.2,
        handlelength=1.4,
    )

    fig.subplots_adjust(left=0.18, right=0.86, top=0.97, bottom=0.16)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    harmonized, exact_short, total = load_agreement()
    mapping = build_mapping_table()
    build_figure(mapping, harmonized, exact_short, total)


if __name__ == "__main__":
    main()
