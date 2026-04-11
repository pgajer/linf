#!/usr/bin/env python3
"""Build manuscript-facing external-validation assets."""

from __future__ import annotations

import math
import re
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from paper_paths import DCST_VALIDATION_DIR, FIGURES_DIR, TABLES_DIR

REBUILT_ROOT = DCST_VALIDATION_DIR / "absorb"
FROZEN_ROOT = DCST_VALIDATION_DIR / "frozen_agp"

TABLE_TSV = TABLES_DIR / "TABLE_V1_external_validation_summary.tsv"
TABLE_MD = TABLES_DIR / "TABLE_V1_external_validation_summary.md"
FIGURE_PNG = FIGURES_DIR / "FIGURE_4_external_validation_modes.png"

KEY_COMPARISONS = [
    {
        "cohort": "halfvarson_2017",
        "comparison": "IBD_vs_Healthy",
        "display": "Halfvarson 2017\nIBD vs Healthy",
        "short": "Halfvarson IBD",
    },
    {
        "cohort": "halfvarson_2017",
        "comparison": "Crohn_vs_Healthy",
        "display": "Halfvarson 2017\nCrohn vs Healthy",
        "short": "Halfvarson Crohn",
    },
    {
        "cohort": "halfvarson_2017",
        "comparison": "UC_vs_Healthy",
        "display": "Halfvarson 2017\nUC vs Healthy",
        "short": "Halfvarson UC",
    },
    {
        "cohort": "hmp2",
        "comparison": "IBD_vs_Healthy",
        "display": "HMP2 / IBDMDB\nIBD vs Healthy",
        "short": "HMP2 IBD",
    },
    {
        "cohort": "gevers_2014",
        "comparison": "Crohn_vs_Healthy",
        "display": "Gevers 2014\nCrohn vs Healthy",
        "short": "Gevers Crohn",
    },
    {
        "cohort": "prjeb84421",
        "comparison": "OFG_vs_Healthy",
        "display": "PRJEB84421\nOFG vs Healthy",
        "short": "PRJEB84421 OFG",
    },
    {
        "cohort": "jacobs_2023_ibs_250bp",
        "comparison": "IBS_vs_Healthy",
        "display": "Jacobs 2023 IBS\nIBS vs Healthy",
        "short": "Jacobs IBS",
    },
]

COHORT_COLORS = {
    "halfvarson_2017": "#1f77b4",
    "hmp2": "#d62728",
    "gevers_2014": "#9467bd",
    "prjeb84421": "#2ca02c",
    "jacobs_2023_ibs_250bp": "#ff7f0e",
}

LINE_STYLES = {
    "IBD_vs_Healthy": "-",
    "Crohn_vs_Healthy": "--",
    "UC_vs_Healthy": ":",
    "OFG_vs_Healthy": "-.",
    "IBS_vs_Healthy": (0, (3, 1, 1, 1)),
}

MARKERS = {
    "IBD_vs_Healthy": "o",
    "Crohn_vs_Healthy": "s",
    "UC_vs_Healthy": "^",
    "OFG_vs_Healthy": "D",
    "IBS_vs_Healthy": "P",
}


def result_path(root: Path, cohort: str, comparison: str, depth: int) -> Path:
    return root / cohort / f"{cohort}__{comparison}_depth{depth}_results.csv"


def mapping_path(cohort: str, comparison: str) -> Path:
    return FROZEN_ROOT / cohort / f"{cohort}__{comparison}_mapping_summary.csv"


def summary_path(root: Path, cohort: str, comparison: str) -> Path:
    return root / cohort / f"{cohort}__{comparison}_validation_summary.md"


def read_best_signal(root: Path, cohort: str, comparison: str) -> dict[str, object]:
    best_row = None
    best_depth = None
    best_q = None
    for depth in range(1, 5):
      path = result_path(root, cohort, comparison, depth)
      if not path.exists():
          continue
      df = pd.read_csv(path)
      if df.empty:
          continue
      row = df.sort_values(["q_value", "p_value", "DCST"]).iloc[0]
      q_val = float(row["q_value"])
      if best_row is None or q_val < best_q:
          best_row = row
          best_depth = depth
          best_q = q_val

    if best_row is None:
        return {
            "depth": np.nan,
            "state": "none",
            "or": np.nan,
            "q_value": np.nan,
            "p_value": np.nan,
            "mapped_samples": np.nan,
            "total_samples": np.nan,
        }

    mapped_samples = float(best_row["mapped_samples"]) if "mapped_samples" in best_row.index else np.nan
    total_samples = (
        float(best_row["total_samples_in_comparison"])
        if "total_samples_in_comparison" in best_row.index
        else np.nan
    )

    return {
        "depth": best_depth,
        "state": str(best_row["DCST"]),
        "or": float(best_row["OR"]) if pd.notna(best_row["OR"]) else np.nan,
        "q_value": float(best_row["q_value"]),
        "p_value": float(best_row["p_value"]),
        "mapped_samples": mapped_samples,
        "total_samples": total_samples,
    }


def read_mapping_coverage(cohort: str, comparison: str) -> pd.DataFrame:
    df = pd.read_csv(mapping_path(cohort, comparison))
    df["cohort"] = cohort
    df["comparison"] = comparison
    return df


def pretty_state(label: str) -> str:
    if not label or label == "none":
        return "none"
    label = label.replace("__", " / ")
    label = label.replace("[Ruminococcus]_torques_group", "Ruminococcus torques group")
    label = label.replace("[Eubacterium]_eligens_group", "Eubacterium eligens group")
    return label


def pretty_q(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 1e-3:
        return f"{value:.2e}"
    return f"{value:.3f}"


def pretty_or(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if math.isinf(value):
        return "Inf"
    return f"{value:.2f}"


def interpret_row(rebuilt_q: float, frozen_q: float) -> str:
    rebuilt_sig = pd.notna(rebuilt_q) and rebuilt_q < 0.05
    frozen_sig = pd.notna(frozen_q) and frozen_q < 0.05
    if rebuilt_sig and frozen_sig:
        return "reproducible and portable"
    if rebuilt_sig and not frozen_sig:
        return "cohort-rebuilt replication only"
    if not rebuilt_sig and frozen_sig:
        return "portable AGP transfer despite weak rebuilt screen"
    return "directional or null external support"


def build_table() -> tuple[pd.DataFrame, pd.DataFrame]:
    rows = []
    coverage_rows = []
    for spec in KEY_COMPARISONS:
        rebuilt = read_best_signal(REBUILT_ROOT, spec["cohort"], spec["comparison"])
        frozen = read_best_signal(FROZEN_ROOT, spec["cohort"], spec["comparison"])
        coverage = read_mapping_coverage(spec["cohort"], spec["comparison"])
        depth2 = coverage.loc[coverage["depth"] == 2]
        depth2_cov = float(depth2["mapping_rate_total"].iloc[0]) if not depth2.empty else np.nan
        rows.append(
            {
                "cohort": spec["short"],
                "comparison": spec["comparison"],
                "rebuilt_best_depth": rebuilt["depth"],
                "rebuilt_best_state": pretty_state(rebuilt["state"]),
                "rebuilt_best_or": pretty_or(rebuilt["or"]),
                "rebuilt_best_q": pretty_q(rebuilt["q_value"]),
                "frozen_best_depth": frozen["depth"],
                "frozen_best_state": pretty_state(frozen["state"]),
                "frozen_best_or": pretty_or(frozen["or"]),
                "frozen_best_q": pretty_q(frozen["q_value"]),
                "frozen_depth2_coverage_pct": f"{100 * depth2_cov:.1f}" if pd.notna(depth2_cov) else "NA",
                "interpretation": interpret_row(rebuilt["q_value"], frozen["q_value"]),
                "display": spec["display"],
                "cohort_slug": spec["cohort"],
                "rebuilt_q_numeric": rebuilt["q_value"],
                "frozen_q_numeric": frozen["q_value"],
            }
        )
        coverage_rows.append(coverage.assign(display=spec["display"], short=spec["short"]))

    return pd.DataFrame(rows), pd.concat(coverage_rows, ignore_index=True)


def write_table(table: pd.DataFrame) -> None:
    out = table[
        [
            "cohort",
            "comparison",
            "rebuilt_best_depth",
            "rebuilt_best_state",
            "rebuilt_best_or",
            "rebuilt_best_q",
            "frozen_best_depth",
            "frozen_best_state",
            "frozen_best_or",
            "frozen_best_q",
            "frozen_depth2_coverage_pct",
            "interpretation",
        ]
    ].copy()
    out.rename(
        columns={
            "rebuilt_best_depth": "rebuilt_depth",
            "rebuilt_best_state": "rebuilt_best_state",
            "rebuilt_best_or": "rebuilt_best_or",
            "rebuilt_best_q": "rebuilt_best_q",
            "frozen_best_depth": "frozen_depth",
            "frozen_best_state": "frozen_best_state",
            "frozen_best_or": "frozen_best_or",
            "frozen_best_q": "frozen_best_q",
            "frozen_depth2_coverage_pct": "frozen_depth2_coverage_pct",
        },
        inplace=True,
    )
    out.to_csv(TABLE_TSV, sep="\t", index=False)

    headers = list(out.columns)
    lines = [
        "| " + " | ".join(headers) + " |",
        "| " + " | ".join(["---"] * len(headers)) + " |",
    ]
    for _, row in out.iterrows():
        vals = [str(row[h]) for h in headers]
        lines.append("| " + " | ".join(vals) + " |")
    TABLE_MD.write_text(
        "# Table V1: External validation summary\n\n" + "\n".join(lines) + "\n",
        encoding="utf-8",
    )


def q_to_score(value: float) -> float:
    if pd.isna(value):
        return 0.0
    value = max(value, 1e-12)
    return -math.log10(value)


def build_figure(table: pd.DataFrame, coverage: pd.DataFrame) -> None:
    order = table.index.tolist()
    y = np.arange(len(order))[::-1]

    fig = plt.figure(figsize=(13.5, 9))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.15, 1], hspace=0.32)
    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[1, 0])

    rebuilt_x = [q_to_score(table.loc[i, "rebuilt_q_numeric"]) for i in order]
    frozen_x = [q_to_score(table.loc[i, "frozen_q_numeric"]) for i in order]
    labels = [table.loc[i, "display"] for i in order]

    for yi, rx, fx in zip(y, rebuilt_x, frozen_x):
        ax1.plot([rx, fx], [yi, yi], color="#b0b0b0", linewidth=1.5, zorder=1)

    ax1.scatter(rebuilt_x, y, s=55, color="#4c78a8", label="Rebuilt absorb", zorder=3)
    ax1.scatter(frozen_x, y, s=55, color="#f58518", label="Frozen AGP transfer", zorder=3)
    ax1.axvline(-math.log10(0.05), color="black", linestyle="--", linewidth=1)
    ax1.set_yticks(y)
    ax1.set_yticklabels(labels, fontsize=9)
    ax1.set_xlabel(r"Best signal strength ($-\log_{10}(q)$)")
    ax1.set_title("A. Rebuilt-cohort validation is broader than frozen AGP transfer", fontsize=12)
    ax1.grid(axis="x", alpha=0.25)
    ax1.legend(loc="lower right", frameon=False)

    for yi, fx, i in zip(y, frozen_x, order):
        q_label = table.loc[i, "frozen_best_q"]
        ax1.text(fx + 0.08, yi, f"q={q_label}", va="center", fontsize=8, color="#444444")

    for spec in KEY_COMPARISONS:
        sub = coverage[
            (coverage["cohort"] == spec["cohort"]) & (coverage["comparison"] == spec["comparison"])
        ].sort_values("depth")
        ax2.plot(
            sub["depth"],
            100 * sub["mapping_rate_total"],
            marker=MARKERS.get(spec["comparison"], "o"),
            linewidth=2,
            markersize=5,
            label=spec["short"],
            color=COHORT_COLORS[spec["cohort"]],
            linestyle=LINE_STYLES.get(spec["comparison"], "-"),
        )
    ax2.set_xlim(1, 4)
    ax2.set_xticks([1, 2, 3, 4])
    ax2.set_ylim(0, 105)
    ax2.set_ylabel("Frozen-transfer mapping coverage (%)")
    ax2.set_xlabel("Depth")
    ax2.set_title("B. Most cohorts map well at depth 1-2, then diverge at deeper levels", fontsize=12)
    ax2.grid(alpha=0.25)
    ax2.legend(loc="lower left", ncol=2, frameon=False, fontsize=8)

    fig.savefig(FIGURE_PNG, dpi=220, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    table, coverage = build_table()
    write_table(table)
    build_figure(table, coverage)
    print(f"Wrote {TABLE_TSV}")
    print(f"Wrote {TABLE_MD}")
    print(f"Wrote {FIGURE_PNG}")


if __name__ == "__main__":
    main()
