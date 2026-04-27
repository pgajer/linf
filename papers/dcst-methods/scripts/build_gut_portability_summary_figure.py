#!/usr/bin/env python3
"""Build main-text and supplement gut portability figures for the dCST methods paper."""

from __future__ import annotations

import math
import textwrap
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D


ROOT = Path("/Users/pgajer/current_projects/linf")
METHODS_ROOT = ROOT / "papers" / "dcst-methods"
GUT_ROOT = ROOT / "papers" / "gut-dcst"
RUN_DIR = Path(
    "/Users/pgajer/current_projects/gut_microbiome/outputs/dcst_analysis/runs/"
    "2026-04-26-agp-silva-local-qza-absorb-depth4-n0_50_25_25_25"
)

FIGURES_DIR = METHODS_ROOT / "assets" / "figures"
TABLES_DIR = GUT_ROOT / "assets" / "tables"

OUT_MAIN = FIGURES_DIR / "FIGURE_2_gut_portability_overview.png"
OUT_SUPP = FIGURES_DIR / "FIGURE_S1_gut_validation_detail.png"

VALIDATION_TABLE = TABLES_DIR / "TABLE_V1_external_validation_summary.tsv"
SUBJECT_TABLE = TABLES_DIR / "TABLE_S12_external_subject_level_sensitivity.tsv"
HALFVARSON_FOLLOWUP_TABLE = TABLES_DIR / "TABLE_S13D_halfvarson_clinical_followup_subject_sensitivity.tsv"

FOCAL_OUTCOMES = [
    ("IBD", "IBD", "#8c1c13"),
    ("Autoimmune", "Autoimmune disease", "#355070"),
    ("IBS", "IBS", "#6a994e"),
]

VALIDATION_ORDER = [
    "Halfvarson IBD",
    "Halfvarson Crohn",
    "Halfvarson UC",
    "HMP2 IBD",
    "Gevers Crohn",
]

SUBJECT_MAP = {
    "Halfvarson IBD": ("Halfvarson 2017", "IBD vs healthy"),
    "Halfvarson Crohn": ("Halfvarson 2017", "Crohn vs healthy"),
    "Halfvarson UC": ("Halfvarson 2017", "UC vs healthy"),
    "HMP2 IBD": ("HMP2 / IBDMDB", "IBD vs healthy"),
    "Gevers Crohn": ("Gevers 2014", "Crohn vs healthy"),
}

MODE_COLORS = {
    "rebuilt_cohort": "#4c78a8",
    "agp_label_transfer": "#f58518",
}

MODE_NAMES = {
    "rebuilt_cohort": "Rebuilt",
    "agp_label_transfer": "AGP transfer",
}


def q_to_score(value: float) -> float:
    if pd.isna(value):
        return 0.0
    return -math.log10(max(float(value), 1e-12))


def fmt_q(value: float) -> str:
    if pd.isna(value):
        return "NA"
    value = float(value)
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


def clean_label(value: str) -> str:
    return str(value).replace("__", " / ")


def load_inputs() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    omnibus = pd.read_csv(RUN_DIR / "omnibus_by_depth.tsv", sep="\t")
    omnibus = omnibus[omnibus["depth"].between(1, 4)].copy()

    validation = pd.read_csv(VALIDATION_TABLE, sep="\t")
    for col in ("rebuilt_best_q", "frozen_best_q"):
        validation[col] = pd.to_numeric(validation[col].replace("NA", np.nan), errors="coerce")

    subject = pd.read_csv(SUBJECT_TABLE, sep="\t")
    for col in ("sample_q", "first_subject_q", "bootstrap_consensus_rate"):
        subject[col] = pd.to_numeric(subject[col], errors="coerce")

    followup = pd.read_csv(HALFVARSON_FOLLOWUP_TABLE, sep="\t")
    for col in ("sample_q", "first_subject_q", "bootstrap_consensus_rate"):
        followup[col] = pd.to_numeric(followup[col], errors="coerce")

    return omnibus, validation, subject, followup


def subject_row(subject: pd.DataFrame, cohort_key: str, mode: str) -> pd.Series:
    cohort, comparison = SUBJECT_MAP[cohort_key]
    sub = subject[
        (subject["cohort"] == cohort)
        & (subject["comparison"] == comparison)
        & (subject["mode"] == mode)
    ]
    if sub.empty:
        raise KeyError(f"Missing subject row for {cohort_key} / {mode}")
    return sub.iloc[0]


def plot_agp_screen(ax: plt.Axes, omnibus: pd.DataFrame) -> None:
    for key, display, color in FOCAL_OUTCOMES:
        sub = omnibus[omnibus["Condition"] == key].sort_values("depth")
        ax.plot(
            sub["depth"],
            sub["cramers_v"],
            marker="o",
            linewidth=2.4,
            ms=6,
            color=color,
            label=display,
        )
    ax.set_xticks([1, 2, 3, 4])
    ax.set_xlabel("dCST hierarchy depth")
    ax.set_ylabel("Cramer's V")
    ax.set_ylim(0.05, 0.38)
    ax.set_title("A. AGP condition screen", loc="left", fontweight="bold")
    ax.grid(axis="y", alpha=0.22, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper left")


def plot_transfer_matrix(ax: plt.Axes, validation: pd.DataFrame, subject: pd.DataFrame) -> None:
    validation = validation.set_index("cohort").loc[VALIDATION_ORDER].reset_index()
    columns = [
        ("rebuilt_cohort", "sample", "Rebuilt\nsample"),
        ("rebuilt_cohort", "first", "Rebuilt\nfirst subject"),
        ("agp_label_transfer", "sample", "Transfer\nsample"),
        ("agp_label_transfer", "first", "Transfer\nfirst subject"),
    ]
    y_positions = np.arange(len(VALIDATION_ORDER))[::-1]

    for row_idx, row in validation.iterrows():
        y = y_positions[row_idx]
        for col_idx, (mode, level, _) in enumerate(columns):
            if mode == "rebuilt_cohort":
                sample_q = row["rebuilt_best_q"]
            else:
                sample_q = row["frozen_best_q"]
            diag = subject_row(subject, row["cohort"], mode)
            q_value = sample_q if level == "sample" else diag["first_subject_q"]
            score = q_to_score(q_value)
            size = 80 + min(score, 10) * 34
            facecolor = MODE_COLORS[mode] if q_value < 0.05 else "white"
            ax.scatter(
                col_idx,
                y,
                s=size,
                facecolors=facecolor,
                edgecolors=MODE_COLORS[mode],
                linewidth=1.8,
                zorder=3,
            )

    ax.set_yticks(y_positions)
    ax.set_yticklabels(VALIDATION_ORDER)
    ax.set_xticks(np.arange(len(columns)))
    ax.set_xticklabels([col[2] for col in columns])
    ax.set_xlim(-0.6, len(columns) - 0.4)
    ax.set_ylim(-0.7, len(VALIDATION_ORDER) - 0.3)
    ax.set_title("B. External IBD portability summary", loc="left", fontweight="bold")
    ax.grid(axis="x", alpha=0.18, linewidth=0.6)
    ax.grid(axis="y", alpha=0.18, linewidth=0.6)
    ax.tick_params(axis="x", labelsize=9.2)
    ax.tick_params(axis="y", labelsize=9.8)
    for spine in ax.spines.values():
        spine.set_visible(False)

    ax.text(
        0.0,
        -0.17,
        "Filled circles indicate q < 0.05; open circles indicate q >= 0.05.",
        transform=ax.transAxes,
        fontsize=8.6,
        color="#444444",
    )


def summarize_halfvarson_followup(followup: pd.DataFrame) -> pd.DataFrame:
    rows = []
    display = {
        "uc_extent": "UC extent",
        "cd_location": "Crohn location",
        "cd_calprotectin": "Crohn calprotectin",
    }
    for analysis, analysis_display in display.items():
        for mode in ("Rebuilt", "AGP transfer"):
            sub = followup[(followup["analysis"] == analysis) & (followup["mode"] == mode)]
            if sub.empty:
                continue
            best = sub.loc[sub["sample_q"].idxmin()]
            rows.append(
                {
                    "analysis": analysis_display,
                    "mode": mode,
                    "sample_q": float(best["sample_q"]),
                    "first_q": float(best["first_subject_q"]),
                    "boot": 100 * float(best["bootstrap_consensus_rate"]),
                    "depth": int(best["depth"]),
                    "label": clean_label(best["top_label"]),
                }
            )
    return pd.DataFrame(rows)


def plot_halfvarson_matrix(ax: plt.Axes, followup: pd.DataFrame) -> None:
    summary = summarize_halfvarson_followup(followup)
    analyses = ["UC extent", "Crohn location", "Crohn calprotectin"]
    modes = ["Rebuilt", "AGP transfer"]
    matrix = np.full((len(analyses), len(modes)), np.nan)
    first_q = np.full_like(matrix, np.nan)
    for _, row in summary.iterrows():
        i = analyses.index(row["analysis"])
        j = modes.index(row["mode"])
        matrix[i, j] = row["boot"]
        first_q[i, j] = row["first_q"]

    im = ax.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=60, aspect="auto")
    for i in range(len(analyses)):
        for j in range(len(modes)):
            if np.isfinite(first_q[i, j]) and first_q[i, j] < 0.05:
                ax.scatter(j, i, marker="*", s=170, color="#202020", zorder=3)
    ax.set_yticks(np.arange(len(analyses)))
    ax.set_yticklabels(analyses)
    ax.set_xticks(np.arange(len(modes)))
    ax.set_xticklabels(modes)
    ax.set_title("C. Halfvarson follow-up", loc="left", fontweight="bold")
    ax.tick_params(axis="x", labelsize=9.4)
    ax.tick_params(axis="y", labelsize=9.8)
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(modes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(analyses), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("Bootstrap recurrence (%)")
    ax.text(0.0, -0.17, "Star indicates first-subject q < 0.05.", transform=ax.transAxes, fontsize=8.6, color="#444444")


def build_main_figure(omnibus: pd.DataFrame, validation: pd.DataFrame, subject: pd.DataFrame, followup: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.size": 10.5,
            "axes.titlesize": 11.5,
            "axes.labelsize": 10.5,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )
    fig = plt.figure(figsize=(12.0, 8.2))
    gs = fig.add_gridspec(2, 2, height_ratios=[0.82, 1.25], width_ratios=[1.04, 0.96], hspace=0.62, wspace=0.48)
    ax_a = fig.add_subplot(gs[0, :])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[1, 1])

    plot_agp_screen(ax_a, omnibus)
    plot_transfer_matrix(ax_b, validation, subject)
    plot_halfvarson_matrix(ax_c, followup)

    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_MAIN, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def plot_detailed_validation(ax: plt.Axes, validation: pd.DataFrame, subject: pd.DataFrame) -> None:
    validation = validation.set_index("cohort").loc[VALIDATION_ORDER].reset_index()
    y = np.arange(len(validation))[::-1]
    offsets = {"rebuilt_cohort": 0.16, "agp_label_transfer": -0.16}
    for idx, row in validation.iterrows():
        yy = y[idx]
        for mode, q_col, label_col in (
            ("rebuilt_cohort", "rebuilt_best_q", "rebuilt_best_state"),
            ("agp_label_transfer", "frozen_best_q", "frozen_best_state"),
        ):
            diag = subject_row(subject, row["cohort"], mode)
            y_mode = yy + offsets[mode]
            x_sample = q_to_score(row[q_col])
            x_first = q_to_score(diag["first_subject_q"])
            color = MODE_COLORS[mode]
            ax.plot([x_sample, x_first], [y_mode, y_mode], color=color, linewidth=1.8, alpha=0.65)
            ax.scatter(x_sample, y_mode, s=65 + 200 * float(diag["bootstrap_consensus_rate"]), color=color, edgecolor="white", linewidth=0.8, zorder=3)
            ax.scatter(x_first, y_mode, s=55, facecolors="none", edgecolors=color, linewidth=1.5, zorder=4)
            ax.text(
                max(x_sample, x_first) + 0.18,
                y_mode,
                f"{MODE_NAMES[mode]}: {clean_label(row[label_col])}; q {fmt_q(row[q_col])} -> {fmt_q(diag['first_subject_q'])}",
                va="center",
                ha="left",
                fontsize=8.4,
                color="#3d3d3d",
            )

    ax.axvline(-math.log10(0.05), color="#303030", linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(validation["cohort"])
    ax.set_xlabel(r"Corrected signal strength ($-\log_{10} q$)")
    ax.set_title("A. External IBD validation details: sample-level q -> first-subject q", loc="left", fontweight="bold")
    ax.set_xlim(-0.2, 9.8)
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)


def plot_detailed_followup(ax: plt.Axes, followup: pd.DataFrame) -> None:
    summary = summarize_halfvarson_followup(followup)
    summary["label_text"] = summary.apply(
        lambda row: f"d{row['depth']}; {row['label']}; q {fmt_q(row['sample_q'])} -> {fmt_q(row['first_q'])}; boot {row['boot']:.1f}%",
        axis=1,
    )
    analyses = ["UC extent", "Crohn location", "Crohn calprotectin"]
    modes = ["Rebuilt", "AGP transfer"]
    matrix = np.full((len(analyses), len(modes)), np.nan)
    text = [["" for _ in modes] for _ in analyses]
    for _, row in summary.iterrows():
        i = analyses.index(row["analysis"])
        j = modes.index(row["mode"])
        matrix[i, j] = row["boot"]
        text[i][j] = "\n".join(textwrap.wrap(row["label_text"], width=42, break_long_words=False))

    im = ax.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=60, aspect="auto")
    ax.set_yticks(np.arange(len(analyses)))
    ax.set_yticklabels(analyses)
    ax.set_xticks(np.arange(len(modes)))
    ax.set_xticklabels(modes)
    for i in range(len(analyses)):
        for j in range(len(modes)):
            color = "white" if matrix[i, j] >= 36 else "#222222"
            ax.text(j, i, text[i][j], ha="center", va="center", fontsize=8.2, color=color, linespacing=1.2)
    ax.set_title("B. Halfvarson within-disease follow-up details", loc="left", fontweight="bold")
    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(modes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(analyses), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    cbar = plt.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cbar.set_label("Bootstrap recurrence (%)")


def build_supplement_figure(validation: pd.DataFrame, subject: pd.DataFrame, followup: pd.DataFrame) -> None:
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )
    fig = plt.figure(figsize=(13.2, 9.0))
    gs = fig.add_gridspec(2, 1, height_ratios=[1.1, 0.95], hspace=0.46)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    plot_detailed_validation(ax_a, validation, subject)
    plot_detailed_followup(ax_b, followup)
    FIGURES_DIR.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_SUPP, dpi=300, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def main() -> None:
    omnibus, validation, subject, followup = load_inputs()
    build_main_figure(omnibus, validation, subject, followup)
    build_supplement_figure(validation, subject, followup)


if __name__ == "__main__":
    main()
