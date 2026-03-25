from __future__ import annotations

from pathlib import Path

import matplotlib.colors as mcolors
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import Rectangle


BASE = Path("/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis")
ADJUSTED = BASE / "full_cohort_adjusted_results.csv"
OUT = Path("/Users/pgajer/current_projects/linf/dev/FIGURE_3_adjusted_association_overview.png")


CONDITION_ORDER = [
    "IBD",
    "Autoimmune",
    "IBS",
    "Acid_reflux",
    "Seasonal_allergies",
    "Kidney_disease",
    "CDI",
    "Lung_disease",
    "Migraine",
    "Cardiovascular_disease",
    "Diabetes",
    "Obesity",
]


HEATMAP_LABELS = {
    "IBD": "IBD",
    "Autoimmune": "Autoimmune",
    "IBS": "IBS",
    "Acid_reflux": "Acid\nreflux",
    "Seasonal_allergies": "Seasonal\nallergies",
    "Kidney_disease": "Kidney\ndisease",
    "CDI": "CDI",
    "Lung_disease": "Lung\ndisease",
    "Migraine": "Migraine",
    "Cardiovascular_disease": "CVD",
    "Diabetes": "Diabetes",
    "Obesity": "Obesity",
}


FOREST_LABELS = {
    "IBD": "IBD",
    "Autoimmune": "Autoimmune",
    "IBS": "IBS",
    "Acid_reflux": "Acid reflux",
    "Seasonal_allergies": "Seasonal allergies",
    "Kidney_disease": "Kidney disease",
    "CDI": "CDI",
    "Lung_disease": "Lung disease",
    "Migraine": "Migraine",
    "Cardiovascular_disease": "CVD",
    "Diabetes": "Diabetes",
    "Obesity": "Obesity",
}


def pretty_taxon(x: str) -> str:
    return x.replace("d__", "")


def format_or(x: float) -> str:
    if x >= 10:
        return f"{x:.1f}"
    if x >= 1:
        return f"{x:.2f}"
    return f"{x:.2f}"


def condition_label(cond: str, mode: str = "heatmap") -> str:
    mapping = HEATMAP_LABELS if mode == "heatmap" else FOREST_LABELS
    return mapping.get(cond, cond.replace("_", " "))


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, list[str], pd.Series]:
    df = pd.read_csv(ADJUSTED)
    for col in ["adj_q", "adj_OR", "adj_CI_lo", "adj_CI_hi", "n_DCST"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    sig = df.loc[df["adj_q"] < 0.05].copy()
    summary = sig.groupby("DCST_short").agg(
        n_sig=("Condition", "count"),
        min_q=("adj_q", "min"),
        n_dcst=("n_DCST", "max"),
    )
    summary["score"] = -np.log10(summary["min_q"]) + 0.25 * summary["n_sig"]
    row_order = summary.sort_values(["score", "n_dcst"], ascending=[False, False]).index.tolist()

    plot_df = df.loc[df["DCST_short"].isin(row_order)].copy()
    top_forest = sig.sort_values("adj_q").head(18).copy()
    top_forest["label"] = top_forest.apply(
        lambda r: f"{pretty_taxon(r['DCST_short'])} x {condition_label(r['Condition'], mode='forest')}",
        axis=1,
    )

    sig_counts = (
        df.assign(sig=df["adj_q"] < 0.05)
        .groupby("Condition")["sig"]
        .sum()
        .reindex(CONDITION_ORDER)
        .fillna(0)
        .astype(int)
    )
    return plot_df, top_forest, row_order, sig_counts


def build_figure(plot_df: pd.DataFrame, top_forest: pd.DataFrame, row_order: list[str], sig_counts: pd.Series) -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    cmap = mcolors.LinearSegmentedColormap.from_list(
        "or_green",
        ["#1b9e77", "#fdfdfd", "#d95f02"],
    )
    norm = mcolors.TwoSlopeNorm(vmin=-1.7, vcenter=0.0, vmax=4.2)

    pivot_or = (
        plot_df.pivot(index="DCST_short", columns="Condition", values="adj_OR")
        .reindex(index=row_order, columns=CONDITION_ORDER)
    )
    pivot_q = (
        plot_df.pivot(index="DCST_short", columns="Condition", values="adj_q")
        .reindex(index=row_order, columns=CONDITION_ORDER)
    )
    n_dcst = plot_df.groupby("DCST_short")["n_DCST"].max()

    fig = plt.figure(figsize=(11.2, 10.7))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.55], hspace=0.42)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    ax1.set_title("A. A small recurring set of taxa carries the adjusted signal", loc="left", fontweight="bold", pad=12)
    n_rows = len(row_order)
    n_cols = len(CONDITION_ORDER)
    ax1.set_xlim(-0.5, n_cols - 0.5)
    ax1.set_ylim(n_rows - 0.5, -1.6)

    for i, taxon in enumerate(row_order):
        for j, cond in enumerate(CONDITION_ORDER):
            or_val = pivot_or.loc[taxon, cond]
            q_val = pivot_q.loc[taxon, cond]
            if pd.isna(or_val):
                face = "#f0f0f0"
                edge = "#dddddd"
            elif pd.notna(q_val) and q_val < 0.05:
                face = cmap(norm(np.clip(np.log2(or_val), -1.7, 4.2)))
                edge = "white"
            else:
                face = "#fbfbfa"
                edge = "#dfdfdf"
            ax1.add_patch(
                Rectangle((j - 0.5, i - 0.5), 1.0, 1.0, facecolor=face, edgecolor=edge, linewidth=0.9)
            )
            if pd.notna(q_val) and q_val < 0.05:
                txt = format_or(or_val)
                ax1.text(
                    j,
                    i,
                    txt,
                    ha="center",
                    va="center",
                    fontsize=9.0,
                    color="#1f1f1f",
                    fontweight="bold" if q_val < 0.01 else "normal",
                )

    ax1.set_xticks(range(n_cols))
    ax1.set_xticklabels([condition_label(c, mode="heatmap") for c in CONDITION_ORDER])
    ax1.set_yticks(range(n_rows))
    ax1.set_yticklabels([f"{pretty_taxon(t)} ({int(n_dcst[t]):,})" for t in row_order])
    ax1.tick_params(axis="x", bottom=False, top=False, labelsize=10)
    plt.setp(ax1.get_xticklabels(), rotation=35, ha="right", rotation_mode="anchor")
    ax1.tick_params(axis="y", length=0)
    for spine in ax1.spines.values():
        spine.set_visible(False)

    for j, cond in enumerate(CONDITION_ORDER):
        ax1.text(
            j,
            -1.05,
            str(int(sig_counts[cond])),
            ha="center",
            va="center",
            fontsize=10.4,
            color="#5a5a5a",
        )
    ax1.text(
        -0.92,
        -1.05,
        "# sig",
        ha="right",
        va="center",
        fontsize=10.3,
        color="#5a5a5a",
    )
    sm = cm.ScalarMappable(norm=norm, cmap=cmap)
    sm.set_array([])
    cbar = fig.colorbar(sm, ax=ax1, fraction=0.028, pad=0.02)
    cbar.set_label("log2 adjusted OR", rotation=90)
    cbar.ax.tick_params(labelsize=9.5)

    forest = top_forest.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(forest))
    colors = [cmap(norm(np.clip(np.log2(or_), -1.7, 4.2))) for or_ in forest["adj_OR"]]
    ax2.set_title("B. Top adjusted associations mix rare high-OR states with common moderate shifts", loc="left", fontweight="bold", pad=12)
    ax2.axvline(1.0, color="#4a4a4a", linewidth=1.2)

    for yi, (_, row), color in zip(y, forest.iterrows(), colors):
        left = row["adj_OR"] - row["adj_CI_lo"]
        right = row["adj_CI_hi"] - row["adj_OR"]
        ax2.errorbar(
            row["adj_OR"],
            yi,
            xerr=np.array([[left], [right]]),
            fmt="o",
            color=color,
            ecolor=color,
            elinewidth=2.0,
            markersize=7.5,
            markerfacecolor=color,
            markeredgecolor=color,
            capsize=0,
            zorder=3,
        )
        ax2.text(
            min(row["adj_CI_hi"] * 1.08, 35.0),
            yi,
            f"OR {format_or(row['adj_OR'])}",
            ha="left",
            va="center",
            fontsize=9.8,
            color="#333333",
        )

    ax2.set_yticks(y)
    ax2.set_yticklabels(forest["label"])
    ax2.set_xscale("log")
    ax2.set_xlim(0.22, 40)
    ax2.set_xticks([0.25, 0.5, 1, 2, 4, 8, 16, 32])
    ax2.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax2.set_xlabel("Adjusted odds ratio (95% CI)")
    ax2.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax2.tick_params(axis="y", labelsize=10.2)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)
    fig.subplots_adjust(left=0.27, right=0.92, top=0.95, bottom=0.07)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    plot_df, top_forest, row_order, sig_counts = load_data()
    build_figure(plot_df, top_forest, row_order, sig_counts)


if __name__ == "__main__":
    main()
