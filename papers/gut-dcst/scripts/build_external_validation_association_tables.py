#!/usr/bin/env python3
"""Build supplement tables for validation-dataset appendices."""

from __future__ import annotations

import hashlib
import math
import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from paper_paths import CANONICAL_TRANSFER_ROOT, GUT_MICROBIOME_ROOT, TABLES_DIR
from taxon_formatting import ITALIC_TAXA

ABSORB_ROOT = GUT_MICROBIOME_ROOT / "outputs" / "dcst_validation" / "absorb"
FROZEN_ROOT = CANONICAL_TRANSFER_ROOT

OUT_DIR = TABLES_DIR / "validation_dataset_appendix_tables"
SUMMARY_TSV = TABLES_DIR / "TABLE_S4D_external_ibd_cohort_characterization.tsv"
SUMMARY_TEX = TABLES_DIR / "TABLE_S4D_external_ibd_cohort_characterization.tex"
MASTER_TSV = TABLES_DIR / "TABLE_S14_S18_validation_dataset_appendix_master.tsv"
TABLES_TEX = TABLES_DIR / "TABLE_S14_S18_validation_dataset_appendices.tex"

POSTERIOR_DRAWS = 20000
ITALIC_TAXA_SET = set(ITALIC_TAXA)

MODE_CONFIG = {
    "rebuilt": {
        "root": ABSORB_ROOT,
        "label": "Rebuilt-cohort absorb validation",
        "caption_prefix": "Rebuilt-cohort absorb dCST associations",
        "all_labels_phrase": "All named absorb dCST labels at this depth are shown.",
        "case_control_phrase": "Case/control totals were",
        "suffix_map": {1: "A", 2: "B", 3: "C", 4: "D"},
    },
    "transfer": {
        "root": FROZEN_ROOT,
        "label": "Direct AGP-derived label transfer",
        "caption_prefix": "Direct AGP-derived absorb dCST label-transfer associations",
        "all_labels_phrase": "All mapped AGP-derived absorb dCST labels at this depth are shown.",
        "case_control_phrase": "Mapped case/control totals at this depth were",
        "suffix_map": {1: "E", 2: "F", 3: "G", 4: "H"},
    },
}


@dataclass(frozen=True)
class ComparisonSpec:
    slug: str
    display: str


@dataclass(frozen=True)
class CohortSpec:
    slug: str
    display: str
    appendix_letter: str
    table_group: str
    cohort_note: str
    selection_note: str
    comparisons: tuple[ComparisonSpec, ...]


SUMMARY_COMPARISONS: dict[str, tuple[str, ...]] = {
    "halfvarson_2017": ("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy"),
    "hmp2": ("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy"),
    "gevers_2014": ("Crohn_vs_Healthy",),
}


APPENDIX_COHORTS: tuple[CohortSpec, ...] = (
    CohortSpec(
        slug="halfvarson_2017",
        display="Halfvarson 2017",
        appendix_letter="A",
        table_group="S14",
        cohort_note=(
            "Longitudinal stool cohort. Primary validation retained all AGP-overlap-filtered "
            "samples; one-sample-per-subject sensitivity retained the first available sample "
            "per participant."
        ),
        selection_note="Sample-level primary analysis with first available sample per participant for sensitivity.",
        comparisons=(
            ComparisonSpec("IBD_vs_Healthy", "IBD vs healthy"),
            ComparisonSpec("Crohn_vs_Healthy", "Crohn vs healthy"),
            ComparisonSpec("UC_vs_Healthy", "UC vs healthy"),
        ),
    ),
    CohortSpec(
        slug="hmp2",
        display="HMP2 / IBDMDB",
        appendix_letter="B",
        table_group="S15",
        cohort_note=(
            "Clinically rich mucosal longitudinal cohort. Primary validation retained all "
            "AGP-overlap-filtered samples; one-sample-per-subject sensitivity retained the "
            "first available sample per participant."
        ),
        selection_note="Sample-level primary analysis with first available sample per participant for sensitivity.",
        comparisons=(
            ComparisonSpec("IBD_vs_Healthy", "IBD vs healthy"),
            ComparisonSpec("Crohn_vs_Healthy", "Crohn vs healthy"),
            ComparisonSpec("UC_vs_Healthy", "UC vs healthy"),
            ComparisonSpec("Crohn_vs_UC", "Crohn vs UC"),
        ),
    ),
    CohortSpec(
        slug="gevers_2014",
        display="Gevers 2014",
        appendix_letter="C",
        table_group="S16",
        cohort_note=(
            "Pediatric stool Crohn cohort with low repeated-subject burden after filtering. "
            "Primary validation retained all AGP-overlap-filtered samples."
        ),
        selection_note="All eligible filtered stool runs retained in the primary analysis; repeated-subject burden is low.",
        comparisons=(
            ComparisonSpec("Crohn_vs_Healthy", "Crohn vs healthy"),
        ),
    ),
    CohortSpec(
        slug="prjeb84421",
        display="PRJEB84421",
        appendix_letter="D",
        table_group="S17",
        cohort_note=(
            "Pediatric stool inflammatory cohort spanning Crohn disease and orofacial "
            "granulomatosis. No repeated subjects were detected after filtering."
        ),
        selection_note="All eligible filtered stool samples retained; no repeated-subject sensitivity was required.",
        comparisons=(
            ComparisonSpec("Inflammatory_vs_Healthy", "Inflammatory vs healthy"),
            ComparisonSpec("OFG_vs_Healthy", "OFG vs healthy"),
            ComparisonSpec("Crohn_vs_Healthy", "Crohn vs healthy"),
            ComparisonSpec("Crohn_vs_OFG", "Crohn vs OFG"),
        ),
    ),
    CohortSpec(
        slug="jacobs_2023_ibs_250bp",
        display="Jacobs 2023 IBS 250bp",
        appendix_letter="E",
        table_group="S18",
        cohort_note=(
            "250bp baseline stool IBS validation subset. No repeated subjects were detected "
            "after filtering."
        ),
        selection_note="All eligible filtered baseline stool samples retained; no repeated-subject sensitivity was required.",
        comparisons=(
            ComparisonSpec("IBS_vs_Healthy", "IBS vs healthy"),
        ),
    ),
)


def latex_escape(text: str) -> str:
    replacements = {
        "\\": r"\textbackslash{}",
        "&": r"\&",
        "%": r"\%",
        "$": r"\$",
        "#": r"\#",
        "_": r"\_",
        "{": r"\{",
        "}": r"\}",
        "~": r"\textasciitilde{}",
        "^": r"\textasciicircum{}",
    }
    return "".join(replacements.get(ch, ch) for ch in str(text))


def fmt_num(value: float, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "NA"
    if abs(value) < 0.001 and value != 0:
        return f"{value:.2e}"
    return f"{value:.{digits}g}"


def fmt_interval(center: float, low: float, high: float) -> str:
    if any(not math.isfinite(v) for v in (center, low, high)):
        return "NA"
    return f"{fmt_num(center)} ({fmt_num(low)}--{fmt_num(high)})"


def stable_seed(*parts: str) -> int:
    joined = "|".join(parts).encode("utf-8")
    digest = hashlib.md5(joined).hexdigest()[:8]
    return int(digest, 16)


def compute_bayes_or(
    *,
    cases_in_label: int,
    controls_in_label: int,
    total_cases: int,
    total_controls: int,
    key: str,
) -> tuple[float, float, float, float]:
    case_fail = max(total_cases - cases_in_label, 0)
    control_fail = max(total_controls - controls_in_label, 0)
    rng = np.random.default_rng(stable_seed(key, str(cases_in_label), str(controls_in_label)))
    p_case = rng.beta(cases_in_label + 0.5, case_fail + 0.5, size=POSTERIOR_DRAWS)
    p_control = rng.beta(controls_in_label + 0.5, control_fail + 0.5, size=POSTERIOR_DRAWS)
    odds_case = p_case / (1 - p_case)
    odds_control = p_control / (1 - p_control)
    or_draws = odds_case / odds_control
    return (
        float(np.median(or_draws)),
        float(np.quantile(or_draws, 0.025)),
        float(np.quantile(or_draws, 0.975)),
        float(np.mean(or_draws > 1)),
    )


def normalize_label(raw: str) -> str:
    text = re.sub(r"[^a-z0-9]+", "_", str(raw).strip().lower()).strip("_")
    if "healthy" in text or text == "hc":
        return "healthy"
    if "crohn" in text or text in {"cd", "icd", "icd_r", "icd_nr"}:
        return "crohn"
    if "ulcerative" in text or text == "uc":
        return "uc"
    if text.startswith("ibs") or "irritable" in text:
        return "ibs"
    if "ofg" in text or "orofacial" in text:
        return "ofg"
    if text.startswith("ibd") or "inflammatory_bowel_disease" in text:
        return "ibd"
    if text.startswith("inflammatory"):
        return "inflammatory"
    return text


def classify_phenotype(raw: str) -> str:
    label = normalize_label(raw)
    if label == "healthy":
        return "Healthy"
    if label == "crohn":
        return "CD"
    if label == "uc":
        return "UC"
    return "Other"


def extract_metadata_label_column(meta: pd.DataFrame) -> str:
    for col in ("final_case_control_label", "case_control_hint", "phenotype_group", "phenotype_hint", "disease_subtype_hint"):
        if col in meta.columns:
            return col
    raise KeyError("No comparison label column found in metadata subset")


def extract_phenotype_column(meta: pd.DataFrame) -> str:
    for col in ("phenotype_group", "phenotype_hint", "disease_subtype_hint", "final_case_control_label"):
        if col in meta.columns:
            return col
    raise KeyError("No phenotype-like column found in metadata subset")


def metadata_path(root: Path, cohort: str, comparison: str) -> Path:
    return root / cohort / f"{cohort}__{comparison}_metadata_subset.csv"


def results_path(root: Path, cohort: str, comparison: str, depth: int) -> Path:
    return root / cohort / f"{cohort}__{comparison}_depth{depth}_results.csv"


def load_metadata_for_mode(root: Path, cohort: str, comparison: str) -> pd.DataFrame:
    meta = pd.read_csv(metadata_path(root, cohort, comparison))
    if "participant_id_hint" not in meta.columns:
        meta["participant_id_hint"] = meta["run_accession"]
    meta["participant_id_hint"] = meta["participant_id_hint"].fillna(meta["run_accession"]).astype(str)
    return meta


def load_metadata(cohort: str, comparison: str) -> pd.DataFrame:
    for root in (FROZEN_ROOT, ABSORB_ROOT):
        path = metadata_path(root, cohort, comparison)
        if path.exists():
            return load_metadata_for_mode(root, cohort, comparison)
    raise FileNotFoundError(f"No metadata subset found for {cohort} {comparison}")


def comparison_control_token(comparison_slug: str) -> str:
    _, right = comparison_slug.split("_vs_")
    return normalize_label(right)


def compute_comparison_totals(meta: pd.DataFrame, comparison_slug: str) -> dict[str, int]:
    label_col = extract_metadata_label_column(meta)
    normalized = meta[label_col].map(normalize_label)
    control_token = comparison_control_token(comparison_slug)
    control_mask = normalized == control_token
    case_mask = ~control_mask
    participant = meta["participant_id_hint"]
    return {
        "total_samples": int(len(meta)),
        "total_cases": int(case_mask.sum()),
        "total_controls": int(control_mask.sum()),
        "total_subjects": int(participant.nunique()),
        "case_subjects": int(participant[case_mask].nunique()),
        "control_subjects": int(participant[control_mask].nunique()),
        "repeated_subjects": int((participant.value_counts() > 1).sum()),
    }


def phenotype_counts(meta: pd.DataFrame) -> tuple[str, str]:
    phenotype_col = extract_phenotype_column(meta)
    broad = meta[phenotype_col].map(classify_phenotype)
    sample_counts = {key: int(value) for key, value in broad.value_counts().to_dict().items()}
    subject_counts = {
        key: int(value)
        for key, value in meta.groupby(broad)["participant_id_hint"].nunique().to_dict().items()
    }
    sample_text = f"CD {sample_counts.get('CD', 0)} / UC {sample_counts.get('UC', 0)} / healthy {sample_counts.get('Healthy', 0)}"
    subject_text = f"CD {subject_counts.get('CD', 0)} / UC {subject_counts.get('UC', 0)} / healthy {subject_counts.get('Healthy', 0)}"
    return sample_text, subject_text


def repeated_subject_summary(meta: pd.DataFrame) -> tuple[int, int]:
    counts = meta["participant_id_hint"].value_counts()
    return int(meta["participant_id_hint"].nunique()), int((counts > 1).sum())


def build_summary_rows() -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    allowed = SUMMARY_COMPARISONS
    by_slug = {cohort.slug: cohort for cohort in APPENDIX_COHORTS}
    for cohort_slug, comparison_slugs in allowed.items():
        cohort = by_slug[cohort_slug]
        for comp in cohort.comparisons:
            if comp.slug not in comparison_slugs:
                continue
            meta = load_metadata(cohort.slug, comp.slug)
            sample_text, subject_text = phenotype_counts(meta)
            n_subjects, repeated_subjects = repeated_subject_summary(meta)
            rows.append(
                {
                    "cohort": cohort.display,
                    "comparison": comp.display,
                    "n_samples": int(len(meta)),
                    "n_subjects": n_subjects,
                    "repeated_subjects": repeated_subjects,
                    "sample_breakdown": sample_text,
                    "subject_breakdown": subject_text,
                    "selection_strategy": cohort.selection_note,
                    "cohort_note": cohort.cohort_note,
                }
            )
    return rows


def format_component(component: str) -> str:
    escaped = latex_escape(component)
    if component in ITALIC_TAXA_SET:
        return rf"\taxon{{{escaped}}}"
    return escaped


def format_dcst_label(label: str) -> str:
    parts = [part for part in str(label).split("__") if part]
    return r"\dcstsep ".join(format_component(part) for part in parts)


def comparison_note(
    *,
    rows: pd.DataFrame,
    display_map: dict[str, str],
    totals_map: dict[str, dict[str, int]],
    mode: str,
) -> str:
    bits: list[str] = []
    for comp_slug in display_map:
        if comp_slug not in set(rows["comparison"]):
            continue
        if mode == "transfer":
            sub = rows[rows["comparison"] == comp_slug]
            total_cases = int(sub["total_cases"].iloc[0])
            total_controls = int(sub["total_controls"].iloc[0])
        else:
            total_cases = totals_map[comp_slug]["total_cases"]
            total_controls = totals_map[comp_slug]["total_controls"]
        bits.append(f"{display_map[comp_slug]} {total_cases}/{total_controls}")
    return "; ".join(bits)


def build_master_table() -> pd.DataFrame:
    rows: list[pd.DataFrame] = []
    for cohort in APPENDIX_COHORTS:
        for comp in cohort.comparisons:
            meta = load_metadata(cohort.slug, comp.slug)
            totals = compute_comparison_totals(meta, comp.slug)
            for mode, config in MODE_CONFIG.items():
                root = config["root"]
                for depth in range(1, 5):
                    path = results_path(root, cohort.slug, comp.slug, depth)
                    if not path.exists():
                        continue
                    df = pd.read_csv(path)
                    if df.empty:
                        continue
                    df = df.copy()
                    df["mode"] = mode
                    df["cohort"] = cohort.display
                    df["cohort_slug"] = cohort.slug
                    df["comparison"] = comp.slug
                    df["comparison_display"] = comp.display
                    df["depth"] = depth
                    df["n_controls"] = df["n_DCST"] - df["n_cases"]
                    if mode == "transfer":
                        df["total_cases"] = df["mapped_cases"]
                        df["total_controls"] = df["mapped_controls"]
                    else:
                        df["total_cases"] = totals["total_cases"]
                        df["total_controls"] = totals["total_controls"]
                    bayes = df.apply(
                        lambda row: compute_bayes_or(
                            cases_in_label=int(row["n_cases"]),
                            controls_in_label=int(row["n_controls"]),
                            total_cases=int(row["total_cases"]),
                            total_controls=int(row["total_controls"]),
                            key=f"{mode}|{cohort.slug}|{comp.slug}|{depth}|{row['DCST']}",
                        ),
                        axis=1,
                        result_type="expand",
                    )
                    bayes.columns = ["bayes_or", "bayes_cri_lo", "bayes_cri_hi", "bayes_pr_gt1"]
                    df = pd.concat([df, bayes], axis=1)
                    rows.append(df)
    master = pd.concat(rows, ignore_index=True)
    master.sort_values(
        ["cohort", "mode", "depth", "comparison_display", "q_value", "p_value", "DCST"],
        inplace=True,
    )
    return master


def write_summary_tex(summary: pd.DataFrame) -> None:
    lines = [
        r"\begin{center}",
        r"\scriptsize",
        r"\begin{tabularx}{\textwidth}{@{} l l r r Y Y Y @{}}",
        r"\toprule",
        r"Cohort & Comparison & Samples & Subjects & Sample breakdown (CD / UC / healthy) & Subject breakdown (CD / UC / healthy) & Selection strategy \\",
        r"\midrule",
    ]
    for _, row in summary.iterrows():
        lines.append(
            " & ".join(
                [
                    latex_escape(row["cohort"]),
                    latex_escape(row["comparison"]),
                    str(int(row["n_samples"])),
                    f"{int(row['n_subjects'])} ({int(row['repeated_subjects'])} repeated)",
                    latex_escape(row["sample_breakdown"]),
                    latex_escape(row["subject_breakdown"]),
                    latex_escape(row["selection_strategy"]),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\bottomrule", r"\end{tabularx}", r"\end{center}"])
    SUMMARY_TEX.write_text("\n".join(lines) + "\n", encoding="utf-8")


def render_longtable(
    cohort: CohortSpec,
    *,
    mode: str,
    depth: int,
    df: pd.DataFrame,
    totals_map: dict[str, dict[str, int]],
) -> str:
    config = MODE_CONFIG[mode]
    table_id = f"{cohort.table_group}{config['suffix_map'][depth]}"
    display_map = {comp.slug: comp.display for comp in cohort.comparisons}
    note = comparison_note(rows=df, display_map=display_map, totals_map=totals_map, mode=mode)
    lines = [
        rf"\noindent\textbf{{Supplementary Table {table_id}.}} {config['caption_prefix']} in {latex_escape(cohort.display)} at depth {depth}. "
        rf"{config['all_labels_phrase']} "
        rf"{config['case_control_phrase']}: {latex_escape(note)}.",
        "",
        r"\begingroup",
        r"\scriptsize",
        r"\setlength{\tabcolsep}{3pt}",
        r"\setlength{\LTleft}{0pt}",
        r"\setlength{\LTright}{0pt}",
        r"\begin{longtable}{@{} p{0.12\textwidth} p{0.26\textwidth} r r p{0.15\textwidth} p{0.07\textwidth} p{0.07\textwidth} p{0.15\textwidth} @{}}",
        r"\toprule",
        r"Comparison & dCST & Cases & Controls & Frequentist OR (95\% CI) & \(p\) & \(q\) & Bayesian OR (95\% CrI) \\",
        r"\midrule",
        r"\endfirsthead",
        r"\toprule",
        r"Comparison & dCST & Cases & Controls & Frequentist OR (95\% CI) & \(p\) & \(q\) & Bayesian OR (95\% CrI) \\",
        r"\midrule",
        r"\endhead",
        r"\bottomrule",
        r"\endfoot",
    ]
    for _, row in df.iterrows():
        lines.append(
            " & ".join(
                [
                    latex_escape(display_map[str(row["comparison"])]),
                    format_dcst_label(row["DCST"]),
                    str(int(row["n_cases"])),
                    str(int(row["n_controls"])),
                    fmt_interval(float(row["OR"]), float(row["CI_low"]), float(row["CI_high"])),
                    fmt_num(float(row["p_value"])),
                    fmt_num(float(row["q_value"])),
                    fmt_interval(float(row["bayes_or"]), float(row["bayes_cri_lo"]), float(row["bayes_cri_hi"])),
                ]
            )
            + r" \\"
        )
    lines.extend([r"\end{longtable}", r"\endgroup", r"\bigskip", ""])
    return "\n".join(lines)


def cohort_overview_text(cohort: CohortSpec) -> str:
    bits: list[str] = []
    for comp in cohort.comparisons:
        meta = load_metadata(cohort.slug, comp.slug)
        totals = compute_comparison_totals(meta, comp.slug)
        bits.append(
            f"{comp.display}: {totals['total_samples']} samples from {totals['total_subjects']} subjects "
            f"({totals['repeated_subjects']} subjects with repeated samples)"
        )
    return "; ".join(bits) + "."


def write_tables_tex(master: pd.DataFrame) -> None:
    lines = [
        r"\clearpage",
        r"\section*{Supplementary Validation Dataset Appendices}",
        "",
        r"The appendices below organize the external association analyses by validation dataset. "
        r"Each cohort appendix reports rebuilt-cohort absorb dCST associations followed by direct "
        r"AGP-derived label-transfer associations across all analyzed depths. Frequentist columns "
        r"report Fisher odds ratios from the sample-level 2x2 tables, and Bayesian columns report "
        r"the corresponding Jeffreys-prior posterior median odds ratio and 95\% credible interval "
        r"derived from the same case/control counts.",
        "",
    ]
    for cohort in APPENDIX_COHORTS:
        lines.append(rf"\subsection*{{Appendix {cohort.appendix_letter}. {latex_escape(cohort.display)}}}")
        lines.append(latex_escape(cohort.cohort_note))
        lines.append("")
        lines.append(rf"\noindent\textit{{Comparisons analyzed.}} {latex_escape(cohort_overview_text(cohort))}")
        lines.append("")
        totals_map = {
            comp.slug: compute_comparison_totals(load_metadata(cohort.slug, comp.slug), comp.slug)
            for comp in cohort.comparisons
        }
        for mode in ("rebuilt", "transfer"):
            lines.append(rf"\noindent\textit{{{MODE_CONFIG[mode]['label']} tables.}}")
            lines.append("")
            for depth in range(1, 5):
                subset = master[
                    (master["cohort_slug"] == cohort.slug)
                    & (master["mode"] == mode)
                    & (master["depth"] == depth)
                ].copy()
                subset["comparison"] = pd.Categorical(
                    subset["comparison"],
                    categories=[comp.slug for comp in cohort.comparisons],
                    ordered=True,
                )
                subset.sort_values(["comparison", "q_value", "p_value", "DCST"], inplace=True)
                lines.append(render_longtable(cohort, mode=mode, depth=depth, df=subset, totals_map=totals_map))
    TABLES_TEX.write_text("\n".join(lines).rstrip() + "\n", encoding="utf-8")


def write_raw_tsvs(master: pd.DataFrame) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    for cohort in APPENDIX_COHORTS:
        for mode in MODE_CONFIG:
            for depth in range(1, 5):
                subset = master[
                    (master["cohort_slug"] == cohort.slug)
                    & (master["mode"] == mode)
                    & (master["depth"] == depth)
                ].copy()
                subset.to_csv(
                    OUT_DIR / f"{cohort.slug}_{mode}_depth{depth}_associations.tsv",
                    sep="\t",
                    index=False,
                )


def build_tables() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    summary = pd.DataFrame(build_summary_rows())
    summary.to_csv(SUMMARY_TSV, sep="\t", index=False)
    write_summary_tex(summary)

    master = build_master_table()
    master.to_csv(MASTER_TSV, sep="\t", index=False)
    write_raw_tsvs(master)
    write_tables_tex(master)


if __name__ == "__main__":
    build_tables()
