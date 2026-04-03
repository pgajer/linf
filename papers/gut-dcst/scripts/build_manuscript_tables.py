from __future__ import annotations

from pathlib import Path
import math

import pandas as pd

from paper_paths import DCST_ANALYSIS_DIR, TABLES_DIR

ANALYSIS = DCST_ANALYSIS_DIR
OUT = TABLES_DIR

DISCOVERY_ORDER = [
    "IBS",
    "IBD",
    "Autoimmune",
    "Acid_reflux",
    "Cardiovascular_disease",
    "Obesity",
    "Seasonal_allergies",
    "Diabetes",
    "CDI",
    "Lung_disease",
    "Migraine",
    "Kidney_disease",
]

COND_LABEL = {
    "IBS": "IBS",
    "IBD": "IBD",
    "Autoimmune": "Autoimmune",
    "Acid_reflux": "Acid reflux",
    "Cardiovascular_disease": "CVD",
    "Obesity": "Obesity",
    "Seasonal_allergies": "Seasonal allergies",
    "Diabetes": "Diabetes",
    "CDI": "CDI",
    "Lung_disease": "Lung disease",
    "Migraine": "Migraine",
    "Kidney_disease": "Kidney disease",
}


def pretty_dcst(label: str) -> str:
    mapping = {
        "RARE_DOMINANT": "Rare dominant",
        "Unassigned;__;__;__;__;__;__": "Unassigned",
        "d__Bacteria;__;__;__;__;__;__": "Unresolved",
        "d__Escherichia-Shigella": "Escherichia-Shigella",
        "d__Prevotella_9": "Prevotella_9",
        "d__Prevotella_7": "Prevotella_7",
        "d__Pasteurellaceae": "Pasteurellaceae",
        "d__Lachnospiraceae": "Lachnospiraceae",
        "d__Corynebacteriaceae": "Corynebacteriaceae",
    }
    if label in mapping:
        return mapping[label]
    if label.startswith("d__"):
        return label.replace("d__", "")
    return label


def sci_or_fixed(x: float) -> str:
    if pd.isna(x):
        return "--"
    if x == 0:
        return "0"
    if abs(x) < 1e-3:
        return f"{x:.1e}"
    if abs(x) >= 100:
        return f"{x:.0f}"
    return f"{x:.3g}"


def pct(x: float) -> str:
    return f"{x:.1f}"


def or_ci(or_: float, lo: float, hi: float) -> str:
    return f"{or_:.2f} ({lo:.2f}-{hi:.2f})"


def save(df: pd.DataFrame, stem: str, title: str) -> None:
    tsv_path = OUT / f"{stem}.tsv"
    md_path = OUT / f"{stem}.md"
    df.to_csv(tsv_path, sep="\t", index=False)
    md_path.write_text(f"# {title}\n\n" + to_markdown(df) + "\n")


def to_markdown(df: pd.DataFrame) -> str:
    cols = [str(c) for c in df.columns]
    header = "| " + " | ".join(cols) + " |"
    sep = "| " + " | ".join(["---"] * len(cols)) + " |"
    body = []
    for row in df.itertuples(index=False, name=None):
        body.append("| " + " | ".join(str(x) for x in row) + " |")
    return "\n".join([header, sep, *body])


def build_table_1() -> None:
    assign = pd.read_csv(ANALYSIS / "full_cohort_dcst_assignments.csv", usecols=["dcst_depth1_short"])
    counts = assign["dcst_depth1_short"].value_counts().rename_axis("dcst").reset_index(name="n")
    total = int(counts["n"].sum())
    top = counts.head(12).copy()
    top["DCST"] = top["dcst"].map(pretty_dcst)
    top["Share (%)"] = top["n"] / total * 100
    top["Cumulative (%)"] = top["Share (%)"].cumsum()
    top.insert(0, "Rank", range(1, len(top) + 1))
    out = top[["Rank", "DCST", "n", "Share (%)", "Cumulative (%)"]].copy()
    out["Share (%)"] = out["Share (%)"].map(pct)
    out["Cumulative (%)"] = out["Cumulative (%)"].map(pct)
    save(out, "TABLE_1_top_depth1_dcsts", "Table 1: Top depth-1 DCSTs in the full AGP cohort")


def build_table_2() -> None:
    adj = pd.read_csv(ANALYSIS / "full_cohort_adjusted_results.csv")
    clean = pd.read_csv(ANALYSIS / "sensitivity_clean_adjusted_results.csv")
    clean_sig = {(r.DCST_short, r.Condition) for _, r in clean[clean["adj_q"] < 0.05].iterrows()}
    rows = []
    sig = adj[adj["adj_q"] < 0.05]
    for cond in DISCOVERY_ORDER:
        sub = sig[sig["Condition"] == cond].sort_values("adj_q")
        if sub.empty:
            continue
        r = sub.iloc[0]
        direction = "Enriched" if r["adj_OR"] > 1 else "Depleted"
        rows.append(
            {
                "Phenotype": COND_LABEL[cond],
                "Direction": direction,
                "DCST": pretty_dcst(r["DCST_short"]),
                "n_DCST": int(r["n_DCST"]),
                "Adjusted OR (95% CI)": or_ci(r["adj_OR"], r["adj_CI_lo"], r["adj_CI_hi"]),
                "adj q": sci_or_fixed(r["adj_q"]),
                "Retained in clean": "Yes" if (r["DCST_short"], r["Condition"]) in clean_sig else "No",
            }
        )
    out = pd.DataFrame(rows)
    save(out, "TABLE_2_headline_adjusted_associations", "Table 2: Headline adjusted depth-1 associations by phenotype")


def summarize_signal(dcst: str, fisher_or: float, adj_or: float | float("nan")) -> str:
    if pd.isna(adj_or):
        return f"{pretty_dcst(dcst)} ({fisher_or:.2f} -> not estimable)"
    return f"{pretty_dcst(dcst)} ({fisher_or:.2f} -> {adj_or:.2f})"


def build_table_3() -> None:
    uni = pd.read_csv(ANALYSIS / "full_cohort_univariate_results.csv")
    adj = pd.read_csv(ANALYSIS / "full_cohort_adjusted_results.csv")
    clean = pd.read_csv(ANALYSIS / "sensitivity_clean_adjusted_results.csv")
    comp = pd.read_csv(ANALYSIS / "univariate_vs_adjusted_comparison.csv")

    rows = []
    for cond in DISCOVERY_ORDER:
        fisher_count = int((uni.loc[uni["Condition"] == cond, "q_value"] < 0.05).sum())
        adj_count = int((adj.loc[adj["Condition"] == cond, "adj_q"] < 0.05).sum())
        clean_count = int((clean.loc[clean["Condition"] == cond, "adj_q"] < 0.05).sum())

        retained = comp[(comp["Condition"] == cond) & (comp["adj_q"] < 0.05)].sort_values("adj_q")
        retained_text = "--"
        if not retained.empty:
            r = retained.iloc[0]
            retained_text = summarize_signal(r["DCST_short"], r["Fisher_OR"], r["adj_OR"])

        atten = comp[
            (comp["Condition"] == cond)
            & (comp["Fisher_q"] < 0.05)
            & ((comp["adj_q"] >= 0.05) | comp["adj_q"].isna())
        ].copy()
        atten = atten[
            ~atten["DCST_short"].str.contains("Unassigned|Eukaryota|Chloroplast", na=False)
        ].sort_values("Fisher_q")
        atten_text = "--"
        if not atten.empty:
            a = atten.iloc[0]
            atten_text = summarize_signal(a["DCST_short"], a["Fisher_OR"], a["adj_OR"])

        rows.append(
            {
                "Phenotype": COND_LABEL[cond],
                "Fisher q<0.05": fisher_count,
                "Adjusted q<0.05": adj_count,
                "Clean q<0.05": clean_count,
                "Representative retained": retained_text,
                "Representative Fisher-only": atten_text,
            }
        )

    out = pd.DataFrame(rows)
    save(out, "TABLE_3_univariate_vs_adjusted_summary", "Table 3: Disease-wise comparison of univariate, adjusted, and clean results")


def build_table_4() -> None:
    summary = pd.DataFrame(
        [
            {
                "Metric": "Overlapping samples",
                "Value": "30,290",
                "Interpretation": "Common sample set available for direct SILVA-GG2 comparison.",
            },
            {
                "Metric": "Harmonized dominant-genus agreement",
                "Value": "27,739 / 30,290 (91.6%)",
                "Interpretation": "Broad depth-1 landscape is stable after genus-level harmonization.",
            },
            {
                "Metric": "Exact short-label equality",
                "Value": "4,566 / 30,290 (15.1%)",
                "Interpretation": "Most disagreements are relabelings or taxonomic splits rather than ecological contradictions.",
            },
            {
                "Metric": "Named depth-1 states",
                "Value": "SILVA: 40 + rare; GG2: 59 + rare",
                "Interpretation": "GG2 resolves some large SILVA states into more named groups.",
            },
        ]
    )
    remap = pd.DataFrame(
        [
            {
                "SILVA state": "Akkermansia",
                "Main GG2 destination(s)": "Akkermansia",
                "Row share": "99.1%",
                "Reading": "Near one-to-one mapping.",
            },
            {
                "SILVA state": "Faecalibacterium",
                "Main GG2 destination(s)": "Faecalibacterium",
                "Row share": "81.7%",
                "Reading": "Mostly preserved with modest spillover.",
            },
            {
                "SILVA state": "Prevotella_9",
                "Main GG2 destination(s)": "Prevotella",
                "Row share": "90.5%",
                "Reading": "High genus-level concordance despite label change.",
            },
            {
                "SILVA state": "Bacteroides",
                "Main GG2 destination(s)": "Phocaeicola_A; Bacteroides_H_857956; unresolved",
                "Row share": "46.6%; 17.4%; 8.4%",
                "Reading": "Large SILVA state is split across multiple GG2 labels.",
            },
            {
                "SILVA state": "Escherichia-Shigella",
                "Main GG2 destination(s)": "Unresolved",
                "Row share": "99.7%",
                "Reading": "GG2 leaves this dominant state largely unresolved.",
            },
            {
                "SILVA state": "Enterobacteriaceae",
                "Main GG2 destination(s)": "Enterobacteriaceae_A_725029; Enterobacterales_737866",
                "Row share": "92.6%; 4.3%",
                "Reading": "Mostly preserved at family/order level with minor split.",
            },
        ]
    )

    metric_tsv = OUT / "TABLE_4_taxonomy_concordance_summary_A.tsv"
    remap_tsv = OUT / "TABLE_4_taxonomy_concordance_summary_B.tsv"
    md = OUT / "TABLE_4_taxonomy_concordance_summary.md"
    summary.to_csv(metric_tsv, sep="\t", index=False)
    remap.to_csv(remap_tsv, sep="\t", index=False)
    md.write_text(
        "# Table 4: SILVA versus GG2 taxonomy concordance summary\n\n"
        "## Table 4A: Agreement metrics\n\n"
        + to_markdown(summary)
        + "\n\n## Table 4B: Representative remappings\n\n"
        + to_markdown(remap)
        + "\n"
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    build_table_1()
    build_table_2()
    build_table_3()
    build_table_4()


if __name__ == "__main__":
    main()
