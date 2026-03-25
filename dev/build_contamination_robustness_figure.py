from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path("/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis")
FULL = BASE / "full_cohort_adjusted_results.csv"
CLEAN = BASE / "sensitivity_clean_adjusted_results.csv"
OUT = Path("/Users/pgajer/current_projects/linf/dev/FIGURE_4_contamination_robustness.png")


def pretty_taxon(x: str) -> str:
    return x.replace("d__", "")


def pretty_condition(x: str) -> str:
    return x.replace("_", " ")


def load_status() -> pd.DataFrame:
    a = pd.read_csv(FULL)
    b = pd.read_csv(CLEAN)
    key = ["DCST_short", "Condition"]
    m = a.merge(b, on=key, how="outer", suffixes=("_full", "_clean"), indicator=True)
    for c in ["adj_q_full", "adj_q_clean"]:
        m[c] = pd.to_numeric(m[c], errors="coerce")
    m["sig_full"] = m["adj_q_full"] < 0.05
    m["sig_clean"] = m["adj_q_clean"] < 0.05
    m["status"] = m.apply(
        lambda r: "kept"
        if r["sig_full"] and r["sig_clean"]
        else ("lost" if r["sig_full"] else ("gained" if r["sig_clean"] else "ns")),
        axis=1,
    )
    m = m[m["status"] != "ns"].copy()
    m["taxon"] = m["DCST_short"].map(pretty_taxon)
    m["condition_pretty"] = m["Condition"].map(pretty_condition)
    return m


def build_figure(df: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    colors = {"kept": "#5abf90", "lost": "#f4a261", "gained": "#4f7cac"}
    fig = plt.figure(figsize=(10.6, 11.6))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.9, 1.2, 1.0], hspace=0.42)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])
    ax3 = fig.add_subplot(gs[2, 0])

    counts = df["status"].value_counts().reindex(["kept", "lost", "gained"]).fillna(0).astype(int)
    full_sig = int((counts["kept"] + counts["lost"]))
    clean_sig = int((counts["kept"] + counts["gained"]))

    ax1.set_title("A. Overall effect of contaminant-aware filtering", loc="left", fontweight="bold")
    ax1.barh([1], [full_sig], color="#d9d9d9", edgecolor="#404040", height=0.35)
    ax1.barh([0], [clean_sig], color="#b8c4ce", edgecolor="#404040", height=0.35)
    ax1.text(full_sig + 0.4, 1, f"{full_sig} significant before filter", va="center", ha="left", fontsize=11)
    ax1.text(clean_sig + 0.4, 0, f"{clean_sig} significant after filter", va="center", ha="left", fontsize=11)
    ax1.text(0.02, 0.92, f"{counts['kept']} kept", transform=ax1.transAxes, color=colors["kept"], fontsize=12, fontweight="bold")
    ax1.text(0.02, 0.82, f"{counts['lost']} lost", transform=ax1.transAxes, color=colors["lost"], fontsize=12, fontweight="bold")
    ax1.text(0.02, 0.72, f"{counts['gained']} gained", transform=ax1.transAxes, color=colors["gained"], fontsize=12, fontweight="bold")
    ax1.set_xlim(0, max(full_sig, clean_sig) + 6)
    ax1.set_ylim(-0.6, 1.6)
    ax1.set_yticks([])
    ax1.set_xlabel("Significant adjusted associations")
    ax1.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.spines["left"].set_visible(False)

    cond = (
        df.groupby(["condition_pretty", "status"])
        .size()
        .unstack(fill_value=0)
        .reindex(columns=["kept", "lost", "gained"], fill_value=0)
        .sort_values(["kept", "lost", "gained"], ascending=False)
    )
    y = np.arange(len(cond))
    left = np.zeros(len(cond))
    ax2.set_title("B. Which disease sections are most affected?", loc="left", fontweight="bold")
    for status in ["kept", "lost", "gained"]:
        vals = cond[status].to_numpy()
        ax2.barh(y, vals, left=left, color=colors[status], edgecolor="white", label=status.capitalize())
        left += vals
    ax2.set_yticks(y)
    ax2.set_yticklabels(cond.index)
    ax2.tick_params(axis="y", labelsize=10.5)
    ax2.invert_yaxis()
    ax2.set_xlabel("Number of significant associations")
    ax2.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    ax2.legend(frameon=False, ncol=3, loc="upper right")

    ax3.axis("off")
    ax3.set_title("C. Representative signals", loc="left", fontweight="bold", pad=8)
    examples = {
        "Robust after filtering": [
            "Morganella × IBD",
            "Bacteroides × IBD",
            "Prevotella_9 × Seasonal allergies",
            "Bacteroides × Autoimmune",
        ],
        "Lost after filtering": [
            "Streptococcus × Autoimmune",
            "Streptococcus × Acid reflux",
            "Staphylococcus × IBS",
            "Bacteroides × Migraine",
        ],
        "New after filtering": [
            "Prevotella_9 × IBS",
            "Escherichia-Shigella × IBD",
            "Prevotella_9 × CDI",
        ],
    }
    box_specs = [
        (0.02, 0.12, 0.30, 0.74, "Robust after filtering", colors["kept"]),
        (0.35, 0.12, 0.30, 0.74, "Lost after filtering", colors["lost"]),
        (0.68, 0.12, 0.30, 0.74, "New after filtering", colors["gained"]),
    ]
    for x0, y0, w, h, title, color in box_specs:
        rect = plt.Rectangle((x0, y0), w, h, transform=ax3.transAxes, facecolor="#fcfcfb", edgecolor=color, linewidth=1.6)
        ax3.add_patch(rect)
        ax3.text(x0 + 0.02, y0 + h - 0.08, title, transform=ax3.transAxes, fontsize=12, fontweight="bold", color=color)
        items = examples[title]
        for i, item in enumerate(items):
            ax3.text(x0 + 0.03, y0 + h - 0.18 - i * 0.16, f"• {item}", transform=ax3.transAxes, fontsize=11, color="#333333")

    fig.subplots_adjust(left=0.11, right=0.98, top=0.97, bottom=0.07)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    df = load_status()
    build_figure(df)


if __name__ == "__main__":
    main()
