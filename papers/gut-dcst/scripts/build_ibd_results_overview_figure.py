#!/usr/bin/env python3
"""Build the IBD overview figure for the gut-dCST paper."""

from __future__ import annotations

import math

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.lines import Line2D

from paper_paths import CANONICAL_AGP_RUN_DIR, FIGURES_DIR, TABLES_DIR
from taxon_formatting import italicize_taxa_mpl

RUN_DIR = CANONICAL_AGP_RUN_DIR
FOLLOWUP_DIR = RUN_DIR / "label_followup" / "ibd"
VALIDATION_TABLE = TABLES_DIR / "TABLE_V1_external_validation_summary.tsv"
SUBJECT_SENSITIVITY_TABLE = TABLES_DIR / "TABLE_S12_external_subject_level_sensitivity.tsv"
HALFVARSON_FOLLOWUP_TABLE = TABLES_DIR / "TABLE_S13D_halfvarson_clinical_followup_subject_sensitivity.tsv"
OUT_FIG = FIGURES_DIR / "FIGURE_2_ibd_results_overview.png"
SUPP_OUT_FIG = FIGURES_DIR / "FIGURE_S1_agp_ibd_followup_states.png"

FOCAL_OUTCOMES = [
    ("IBD", "IBD", "#8c1c13"),
    ("Autoimmune", "Autoimmune disease", "#355070"),
    ("IBS", "IBS", "#6a994e"),
]

IBD_VALIDATION_ORDER = [
    "Halfvarson IBD",
    "Halfvarson Crohn",
    "Halfvarson UC",
    "HMP2 IBD",
    "Gevers Crohn",
]

VALIDATION_SUBJECT_MAP = {
    "Halfvarson IBD": ("Halfvarson 2017", "IBD vs healthy"),
    "Halfvarson Crohn": ("Halfvarson 2017", "Crohn vs healthy"),
    "Halfvarson UC": ("Halfvarson 2017", "UC vs healthy"),
    "HMP2 IBD": ("HMP2 / IBDMDB", "IBD vs healthy"),
    "Gevers Crohn": ("Gevers 2014", "Crohn vs healthy"),
}

FOLLOWUP_ANALYSIS_ORDER = [
    ("uc_extent", "UC extent"),
    ("cd_location", "Crohn location"),
    ("cd_calprotectin", "Crohn calprotectin"),
]

MODE_STYLES = {
    "rebuilt_cohort": {"color": "#4c78a8", "label": "Rebuilt"},
    "agp_label_transfer": {"color": "#f58518", "label": "AGP transfer"},
}


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


def load_subject_sensitivity() -> pd.DataFrame:
    table = pd.read_csv(SUBJECT_SENSITIVITY_TABLE, sep="\t")
    for col in (
        "sample_q",
        "first_subject_q",
        "bootstrap_consensus_rate",
        "bootstrap_median_q",
    ):
        table[col] = pd.to_numeric(table[col], errors="coerce")
    return table


def load_halfvarson_followup() -> pd.DataFrame:
    table = pd.read_csv(HALFVARSON_FOLLOWUP_TABLE, sep="\t")
    for col in ("sample_q", "first_subject_q", "bootstrap_consensus_rate"):
        table[col] = pd.to_numeric(table[col], errors="coerce")
    return table


def q_to_score(value: float) -> float:
    if pd.isna(value):
        return 0.0
    return -math.log10(max(float(value), 1e-12))


def fmt_q_text(value: float) -> str:
    if pd.isna(value):
        return "NA"
    if value < 0.001:
        return f"{value:.1e}"
    return f"{value:.3f}".rstrip("0").rstrip(".")


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


def panel_b_state_label(label: str) -> str:
    return wrap_state_label(label, max_parts_per_line=2)


def followup_sort_score(row: pd.Series) -> float:
    q_val = row.get("adj_freq_q_global", np.nan)
    if pd.notna(q_val):
        return 100.0 + q_to_score(q_val)
    bayes_or = row.get("bayes_adj_or", np.nan)
    if pd.notna(bayes_or) and bayes_or > 0:
        return abs(math.log(float(bayes_or)))
    return -np.inf


def select_supplement_rows(followup: pd.DataFrame, n_rows: int = 7) -> list[pd.Series]:
    df = followup.copy()
    df["score"] = df.apply(followup_sort_score, axis=1)
    df = df[np.isfinite(df["score"])].copy()
    selected: list[pd.Series] = []
    chosen_labels = set()

    for depth in [1, 2, 3, 4]:
        sub = df[df["depth"] == depth].sort_values(
            ["score", "n_label_all"], ascending=[False, False]
        )
        if sub.empty:
            continue
        row = None
        for _, cand in sub.iterrows():
            if str(cand["label"]) not in chosen_labels:
                row = cand
                break
        if row is None:
            row = sub.iloc[0]
        selected.append(row)
        chosen_labels.add(str(row["label"]))

    extras = df.sort_values(["score", "n_label_all"], ascending=[False, False])
    for _, row in extras.iterrows():
        if str(row["label"]) in chosen_labels:
            continue
        selected.append(row)
        chosen_labels.add(str(row["label"]))
        if len(selected) >= n_rows:
            break

    selected = selected[:n_rows]
    selected.sort(key=lambda row: (int(row["depth"]), -followup_sort_score(row)))
    return selected


def point_size_from_bootstrap(rate: float) -> float:
    if pd.isna(rate):
        rate = 0.0
    return 40.0 + 150.0 * float(rate)


def get_subject_row(subject: pd.DataFrame, cohort_key: str, mode: str) -> pd.Series:
    cohort, comparison = VALIDATION_SUBJECT_MAP[cohort_key]
    sub = subject[
        (subject["cohort"] == cohort)
        & (subject["comparison"] == comparison)
        & (subject["mode"] == mode)
    ]
    if sub.empty:
        raise KeyError(f"Missing subject sensitivity row for {cohort_key} / {mode}")
    return sub.iloc[0]


def summarize_followup_cells(followup: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for analysis_key, analysis_display in FOLLOWUP_ANALYSIS_ORDER:
        for mode_key, style in MODE_STYLES.items():
            sub = followup[
                (followup["analysis"] == analysis_key)
                & (followup["mode"] == style["label"])
            ].copy()
            if sub.empty:
                continue
            best = sub.loc[sub["sample_q"].idxmin()]
            rows.append(
                {
                    "analysis": analysis_display,
                    "mode": style["label"],
                    "depth": int(best["depth"]),
                    "sample_q": float(best["sample_q"]),
                    "first_subject_q": float(best["first_subject_q"]),
                    "bootstrap_rate_pct": 100 * float(best["bootstrap_consensus_rate"]),
                }
            )
    return pd.DataFrame(rows)


def plot_validation_preview(
    ax: plt.Axes, validation: pd.DataFrame, subject: pd.DataFrame
) -> None:
    val = validation[validation["cohort"].isin(IBD_VALIDATION_ORDER)].copy()
    val["sort_order"] = val["cohort"].map(
        {name: idx for idx, name in enumerate(IBD_VALIDATION_ORDER)}
    )
    val.sort_values("sort_order", inplace=True)
    val.reset_index(drop=True, inplace=True)

    y = np.arange(len(val))[::-1]
    offsets = {"rebuilt_cohort": 0.16, "agp_label_transfer": -0.16}
    x_values = []

    for idx, row in val.iterrows():
        yy = y[idx]
        for mode, q_col in (
            ("rebuilt_cohort", "rebuilt_best_q"),
            ("agp_label_transfer", "frozen_best_q"),
        ):
            diag = get_subject_row(subject, row["cohort"], mode)
            x_sample = q_to_score(row[q_col])
            x_first = q_to_score(diag["first_subject_q"])
            x_values.extend([x_sample, x_first])
            y_mode = yy + offsets[mode]
            color = MODE_STYLES[mode]["color"]

            ax.plot(
                [x_sample, x_first],
                [y_mode, y_mode],
                color=color,
                linewidth=1.6,
                alpha=0.6,
                zorder=1,
            )
            ax.scatter(
                [x_sample],
                [y_mode],
                s=point_size_from_bootstrap(diag["bootstrap_consensus_rate"]),
                color=color,
                edgecolor="white",
                linewidth=0.9,
                zorder=3,
            )
            ax.scatter(
                [x_first],
                [y_mode],
                s=44,
                facecolors="none",
                edgecolors=color,
                linewidth=1.3,
                zorder=4,
            )

        label_x = max(q_to_score(row["frozen_best_q"]) + 0.2, 0.22)
        if yy <= 0.25:
            label_y = yy + 0.28
            label_va = "bottom"
        else:
            label_y = yy - 0.27
            label_va = "top"
        ax.text(
            label_x,
            label_y,
            panel_b_state_label(row["frozen_best_state"]),
            ha="left",
            va=label_va,
            fontsize=7.8,
            color="#4b4b4b",
        )

    ax.axvline(-math.log10(0.05), color="#303030", linestyle="--", linewidth=1.0)
    ax.set_yticks(y)
    ax.set_yticklabels(val["cohort"], fontsize=9)
    ax.set_xlabel(r"Corrected signal strength ($-\log_{10} q$)")
    ax.set_title(
        "B. External IBD transfer remains strongest in Halfvarson; open markers show first-subject attenuation and filled-marker area scales with bootstrap support",
        loc="left",
        fontweight="bold",
    )
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.set_xlim(left=-0.35, right=max(max(x_values) + 2.2, 8.0))

    legend_handles = [
        Line2D(
            [0],
            [0],
            color=MODE_STYLES["rebuilt_cohort"]["color"],
            marker="o",
            linewidth=1.6,
            markersize=6,
            label="Rebuilt sample-level",
        ),
        Line2D(
            [0],
            [0],
            color=MODE_STYLES["agp_label_transfer"]["color"],
            marker="o",
            linewidth=1.6,
            markersize=6,
            label="AGP-transfer sample-level",
        ),
        Line2D(
            [0],
            [0],
            color="#555555",
            marker="o",
            markerfacecolor="none",
            linewidth=0,
            markersize=6,
            label="First subject per participant",
        ),
        plt.scatter([], [], s=point_size_from_bootstrap(0.15), color="#7f7f7f", label="Boot 15%"),
        plt.scatter([], [], s=point_size_from_bootstrap(0.30), color="#7f7f7f", label="Boot 30%"),
        plt.scatter([], [], s=point_size_from_bootstrap(0.50), color="#7f7f7f", label="Boot 50%"),
    ]
    ax.legend(handles=legend_handles, frameon=False, loc="lower right", ncol=2, fontsize=8.0)


def plot_halfvarson_followup(ax: plt.Axes, followup: pd.DataFrame) -> None:
    summary = summarize_followup_cells(followup)
    analyses = [display for _, display in FOLLOWUP_ANALYSIS_ORDER]
    modes = [MODE_STYLES["rebuilt_cohort"]["label"], MODE_STYLES["agp_label_transfer"]["label"]]

    matrix = np.full((len(analyses), len(modes)), np.nan)
    text_matrix = [["" for _ in modes] for _ in analyses]

    for _, row in summary.iterrows():
        i = analyses.index(row["analysis"])
        j = modes.index(row["mode"])
        matrix[i, j] = row["bootstrap_rate_pct"]
        text_matrix[i][j] = (
            f"d{row['depth']}\n"
            f"{fmt_q_text(row['sample_q'])} → {fmt_q_text(row['first_subject_q'])}\n"
            f"boot {row['bootstrap_rate_pct']:.1f}%"
        )

    im = ax.imshow(matrix, cmap="YlOrRd", vmin=0, vmax=50, aspect="auto")
    ax.set_xticks(np.arange(len(modes)))
    ax.set_xticklabels(modes)
    ax.set_yticks(np.arange(len(analyses)))
    ax.set_yticklabels(analyses)
    ax.set_title(
        "C. Within Halfvarson, rebuilt-mode Crohn location comes closest to the working bootstrap-stability threshold; cell text reports strongest sample-level $q \\rightarrow$ first-subject $q$",
        loc="left",
        fontweight="bold",
    )
    ax.text(
        0.02,
        1.02,
        "cell = d / sample q \u2192 first-subject q / boot %",
        transform=ax.transAxes,
        fontsize=8.1,
        color="#4f4f4f",
        va="bottom",
    )

    for i in range(len(analyses)):
        for j in range(len(modes)):
            if np.isnan(matrix[i, j]):
                continue
            text_color = "#222222" if matrix[i, j] < 30 else "white"
            ax.text(j, i, text_matrix[i][j], ha="center", va="center", fontsize=8.2, color=text_color)

    cbar = plt.colorbar(im, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("Bootstrap recurrence (%)")

    for spine in ax.spines.values():
        spine.set_visible(False)
    ax.set_xticks(np.arange(-0.5, len(modes), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(analyses), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=1.5)
    ax.tick_params(which="minor", bottom=False, left=False)


def build_main_figure(
    omnibus: pd.DataFrame,
    validation: pd.DataFrame,
    subject: pd.DataFrame,
    followup: pd.DataFrame,
) -> None:
    fig = plt.figure(figsize=(13.4, 11.0))
    gs = fig.add_gridspec(3, 1, height_ratios=[0.9, 1.2, 1.0], hspace=0.58)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[1, 0])
    ax_c = fig.add_subplot(gs[2, 0])

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
    ax_a.set_title(
        "A. In the large AGP training cohort, self-reported IBD is the strongest depth-resolved signal",
        loc="left",
        fontweight="bold",
    )
    ax_a.grid(axis="y", alpha=0.22, linewidth=0.6)
    ax_a.spines["top"].set_visible(False)
    ax_a.spines["right"].set_visible(False)
    ax_a.legend(frameon=False, loc="upper left")

    plot_validation_preview(ax_b, validation, subject)
    plot_halfvarson_followup(ax_c, followup)

    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG, dpi=220, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def build_supplement_figure(followup: pd.DataFrame) -> None:
    fig, ax = plt.subplots(figsize=(11.2, 7.8))
    rows = select_supplement_rows(followup)
    y = np.arange(len(rows))[::-1]
    labels = [
        wrap_state_label(str(row["label"]).replace("__", " / ")) + "\n" + support_label(row, int(row["depth"]))
        for row in rows
    ]
    for idx, row in enumerate(rows):
        yy = y[idx]
        if np.isfinite(row["freq_adj_or"]):
            ax.plot(
                [row["freq_adj_ci_lo"], row["freq_adj_ci_hi"]],
                [yy + 0.12, yy + 0.12],
                color="#355070",
                linewidth=2.0,
            )
            ax.scatter(
                row["freq_adj_or"],
                yy + 0.12,
                color="#355070",
                s=34,
                zorder=3,
                label="Frequentist adjusted OR" if idx == 0 else None,
            )
        ax.plot(
            [row["bayes_adj_ci_lo"], row["bayes_adj_ci_hi"]],
            [yy - 0.12, yy - 0.12],
            color="#bc6c25",
            linewidth=2.0,
        )
        ax.scatter(
            row["bayes_adj_or"],
            yy - 0.12,
            color="#bc6c25",
            s=34,
            zorder=3,
            label="Bayesian adjusted OR" if idx == 0 else None,
        )
    ax.axvline(1.0, color="#6c757d", linestyle="--", linewidth=1.0)
    ax.set_xscale("log")
    ax.set_xlabel("Adjusted odds ratio for self-reported AGP IBD")
    ax.set_yticks(y)
    ax.set_yticklabels(labels)
    ax.tick_params(axis="y", labelsize=8.2)
    ax.grid(axis="x", alpha=0.22, linewidth=0.6)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.legend(frameon=False, loc="upper right", fontsize=8.5)

    SUPP_OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(SUPP_OUT_FIG, dpi=220, bbox_inches="tight", pad_inches=0.05)
    plt.close(fig)


def build_figure() -> None:
    omnibus = load_omnibus()
    followup = load_ibd_followup()
    validation = load_validation()
    subject = load_subject_sensitivity()
    halfvarson_followup = load_halfvarson_followup()

    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.titlesize": 12,
            "axes.labelsize": 10,
            "figure.facecolor": "white",
            "axes.facecolor": "#fcfcfb",
        }
    )
    build_main_figure(omnibus, validation, subject, halfvarson_followup)
    build_supplement_figure(followup)


if __name__ == "__main__":
    build_figure()
