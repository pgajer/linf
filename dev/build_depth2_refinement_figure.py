from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


BASE = Path("/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis")
ASSIGN = BASE / "full_cohort_dcst_assignments.csv"
DEPTH2 = BASE / "depth2_adjusted_results.csv"
DEPTH1 = BASE / "full_cohort_adjusted_results.csv"
OUT = Path("/Users/pgajer/current_projects/linf/dev/FIGURE_6_depth2_refinement_case_study.png")


def parse_last_label(part: str) -> str:
    best = None
    for token in part.split(";"):
        if "__" not in token:
            continue
        _, value = token.split("__", 1)
        if value and value != "_":
            best = value
    return best or "Unknown"


def parse_pair(dcst: str) -> tuple[str, str]:
    parts = dcst.split(";s____")
    if len(parts) == 1:
        parts = dcst.split("____")
    labels = [parse_last_label(part) for part in parts]
    if len(labels) == 1:
        return labels[0], labels[0]
    return labels[0], labels[1]


def color_for_row(adj_or: float, adj_q: float) -> str:
    if adj_q < 0.05 and adj_or > 1:
        return "#d95f02"
    if adj_q < 0.05 and adj_or < 1:
        return "#1b9e77"
    if adj_or > 1.1:
        return "#fdbb84"
    if adj_or < 0.9:
        return "#a6dba0"
    return "#bdbdbd"


def load_data() -> tuple[pd.DataFrame, pd.Series]:
    assign = pd.read_csv(ASSIGN, usecols=["dcst_depth2"])
    parsed = assign["dcst_depth2"].map(parse_pair)
    assign["parent"] = parsed.map(lambda x: x[0])
    assign["child"] = parsed.map(lambda x: x[1])
    counts = assign.loc[assign["parent"] == "Bacteroides", "child"].value_counts()

    res = pd.read_csv(DEPTH2)
    parsed = res["DCST"].map(parse_pair)
    res["parent"] = parsed.map(lambda x: x[0])
    res["child"] = parsed.map(lambda x: x[1])
    for col in ["adj_OR", "adj_CI_lo", "adj_CI_hi", "adj_q"]:
        res[col] = pd.to_numeric(res[col], errors="coerce")
    ibd = res.loc[(res["Condition"] == "IBD") & (res["parent"] == "Bacteroides")].copy()

    total_parent = int(counts.sum())
    top = counts.head(10).rename("count").reset_index().rename(columns={"index": "child"})
    top["share"] = top["count"] / counts.sum()
    top["label"] = top["child"].replace({"Prevotella_9": "Prevotella_9"})
    top = top.merge(
        ibd[["child", "adj_OR", "adj_CI_lo", "adj_CI_hi", "adj_q"]],
        on="child",
        how="left",
        validate="one_to_one",
    )
    top["color"] = [
        color_for_row(or_, q_) for or_, q_ in zip(top["adj_OR"], top["adj_q"])
    ]
    top["display_label"] = top["child"]

    depth1 = pd.read_csv(DEPTH1)
    for col in ["adj_OR", "adj_CI_lo", "adj_CI_hi", "adj_q"]:
        depth1[col] = pd.to_numeric(depth1[col], errors="coerce")
    parent_row = depth1.loc[
        (depth1["Condition"] == "IBD") & (depth1["DCST_short"] == "d__Bacteroides")
    ].iloc[0]

    top.attrs["total_parent"] = total_parent
    return top, parent_row


def build_figure(top: pd.DataFrame, parent_row: pd.Series) -> None:
    plt.rcParams.update(
        {
            "font.size": 11,
            "axes.titlesize": 13,
            "axes.labelsize": 11,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )

    fig = plt.figure(figsize=(11.2, 9.4))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.08, 1.28], hspace=0.42)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    top = top.copy()
    total_parent = int(top.attrs["total_parent"])
    other_count = total_parent - int(top["count"].sum())
    other_share = other_count / total_parent
    plot_a = pd.concat(
        [
            top[["display_label", "count", "share", "color"]],
            pd.DataFrame(
                [{"display_label": "Other 17 children", "count": other_count, "share": other_share, "color": "#d9d9d9"}]
            ),
        ],
        ignore_index=True,
    )
    plot_a = plot_a.iloc[::-1].reset_index(drop=True)
    y = np.arange(len(plot_a))
    ax1.set_title("A. The broad Bacteroides parent contains many common depth-2 children", loc="left", fontweight="bold")
    ax1.barh(y, plot_a["count"], color=plot_a["color"], edgecolor="white", height=0.72)
    ax1.set_yticks(y)
    ax1.set_yticklabels(plot_a["display_label"])
    ax1.set_xlabel("Samples within the depth-1 Bacteroides state")
    ax1.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    for yi, (_, row) in zip(y, plot_a.iterrows()):
        ax1.text(
            row["count"] + 35,
            yi,
            f"{int(row['count']):,} ({row['share']*100:.0f}%)",
            va="center",
            ha="left",
            fontsize=10.5,
        )
    ax1.text(
        0.99,
        0.94,
        f"Top 10 children cover {top['count'].sum()/total_parent:.0%} of all Bacteroides-dominant samples",
        transform=ax1.transAxes,
        fontsize=10.6,
        color="#4a4a4a",
        ha="right",
        va="top",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5, "alpha": 0.9},
    )

    plot_b = top.iloc[::-1].reset_index(drop=True)
    y2 = np.arange(len(plot_b))
    ax2.set_title("B. IBD odds ratios differ materially across those common children", loc="left", fontweight="bold")
    ax2.axvline(1.0, color="#4a4a4a", linewidth=1.2)
    ax2.axvline(parent_row["adj_OR"], color="#7a7a7a", linewidth=1.3, linestyle="--")
    ax2.text(
        parent_row["adj_OR"] * 1.02,
        0.98,
        f"depth-1 parent OR = {parent_row['adj_OR']:.2f}",
        transform=ax2.get_xaxis_transform(),
        ha="left",
        va="top",
        fontsize=10.3,
        color="#5a5a5a",
        bbox={"facecolor": "white", "edgecolor": "none", "pad": 1.5, "alpha": 0.9},
    )
    for yi, (_, row) in zip(y2, plot_b.iterrows()):
        left = row["adj_OR"] - row["adj_CI_lo"]
        right = row["adj_CI_hi"] - row["adj_OR"]
        face = row["color"] if row["adj_q"] < 0.05 else "white"
        ax2.errorbar(
            row["adj_OR"],
            yi,
            xerr=np.array([[left], [right]]),
            fmt="o",
            color=row["color"],
            ecolor=row["color"],
            elinewidth=1.8,
            markersize=7.8,
            markerfacecolor=face,
            markeredgecolor=row["color"],
            markeredgewidth=1.8,
            capsize=0,
            zorder=3,
        )
        q_text = "q<0.05" if row["adj_q"] < 0.05 else "ns"
        ax2.text(
            row["adj_CI_hi"] * 1.08,
            yi,
            f"OR {row['adj_OR']:.2f}, {q_text}",
            va="center",
            ha="left",
            fontsize=10.2,
            color="#333333",
        )
    ax2.set_xscale("log")
    ax2.set_xlim(0.23, 7.2)
    ax2.set_xticks([0.25, 0.5, 1, 2, 4])
    ax2.get_xaxis().set_major_formatter(plt.ScalarFormatter())
    ax2.set_yticks(y2)
    ax2.set_yticklabels(plot_b["display_label"])
    ax2.set_xlabel("Adjusted odds ratio for IBD")
    ax2.grid(axis="x", alpha=0.25, linewidth=0.6)
    ax2.spines["top"].set_visible(False)
    ax2.spines["right"].set_visible(False)

    fig.subplots_adjust(left=0.22, right=0.92, top=0.97, bottom=0.09)
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT, dpi=240, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def main() -> None:
    top, parent_row = load_data()
    build_figure(top, parent_row)


if __name__ == "__main__":
    main()
