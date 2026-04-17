from __future__ import annotations

import gzip
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
PROJECT_ROOT = ROOT.parent
GUT_ROOT = PROJECT_ROOT / "gut_microbiome"

COUNTS_PATH = (
    GUT_ROOT
    / "outputs/prime_species/prime_gut_projects_silva_species_absolute_2026-03-24.csv.gz"
)
METADATA_PATH = GUT_ROOT / "data/prime_gut_project_sample_metadata_2026-03-24.csv.gz"
ASSIGNMENTS_PATH = (
    GUT_ROOT
    / "outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive/agp_absorb_assignments.tsv.gz"
)
OUT = ROOT / "papers/gut-dcst/assets/tables/TABLE_S8_ibd_dominance_lineage_retained_feature_profiles.tsv"

LINEAGES = [
    ("Bacteroides", 1, "Bacteroides"),
    ("Bacteroides__Lachnospiraceae", 2, "Bacteroides/Lachnospiraceae"),
    (
        "Bacteroides__Lachnospiraceae__Alistipes",
        3,
        "Bacteroides/Lachnospiraceae/Alistipes",
    ),
    (
        "Bacteroides__Lachnospiraceae__Alistipes__Faecalibacterium",
        4,
        "Bacteroides/Lachnospiraceae/Alistipes/Faecalibacterium",
    ),
]
COMPONENTS = ["Bacteroides", "Lachnospiraceae", "Alistipes", "Faecalibacterium"]


def extract_rank_value(parts: list[str], prefix: str) -> str:
    for part in parts:
        if part.startswith(prefix):
            value = part[len(prefix) :].strip()
            if value and not set(value) <= {"_"}:
                return value
    return ""


def compact_taxon_label(taxonomy: str) -> str:
    if pd.isna(taxonomy) or not str(taxonomy):
        return "Unassigned"

    parts = str(taxonomy).split(";")
    species = extract_rank_value(parts, "s__")
    genus = extract_rank_value(parts, "g__")
    family = extract_rank_value(parts, "f__")
    order = extract_rank_value(parts, "o__")
    class_name = extract_rank_value(parts, "c__")
    phylum = extract_rank_value(parts, "p__")
    domain = extract_rank_value(parts, "d__")

    if species:
        if genus and not species.startswith(genus):
            return f"{genus} {species}"
        return species
    if genus:
        return genus
    if family:
        return family
    if order:
        return order
    if class_name:
        return class_name
    if phylum:
        return phylum
    if domain:
        return domain
    return "Unassigned"


def make_unique(labels: list[str]) -> list[str]:
    seen: dict[str, int] = {}
    out: list[str] = []
    for label in labels:
        count = seen.get(label, 0)
        out.append(label if count == 0 else f"{label}_{count}")
        seen[label] = count + 1
    return out


def read_csv_gz(path: Path, **kwargs) -> pd.DataFrame:
    with gzip.open(path, "rt") as handle:
        return pd.read_csv(handle, **kwargs)


def main() -> None:
    metadata = read_csv_gz(METADATA_PATH, usecols=["Run", "BioProject"])
    agp_runs = set(metadata.loc[metadata["BioProject"] == "PRJEB11419", "Run"])

    counts = read_csv_gz(COUNTS_PATH)
    counts = counts[counts["Run"].isin(agp_runs)].copy()
    counts = counts.set_index("Run")

    library_size = counts.sum(axis=1)
    counts = counts.loc[library_size >= 1000]
    keep_features = (counts >= 2).sum(axis=0) >= int(0.05 * counts.shape[0])
    counts = counts.loc[:, keep_features]

    labels = make_unique([compact_taxon_label(col) for col in counts.columns])
    counts.columns = labels

    assignments = pd.read_csv(ASSIGNMENTS_PATH, sep="\t")
    assignments = assignments.set_index("Run").loc[counts.index]

    rel = counts.div(counts.sum(axis=1), axis=0)
    rows = []
    for raw_label, depth, display_label in LINEAGES:
        col = f"depth{depth}_absorb"
        mask = assignments[col] == raw_label
        sub_assign = assignments.loc[mask]
        sub_rel = rel.loc[mask, COMPONENTS]
        ibd_rel = sub_rel.loc[sub_assign["IBD"] == 1]

        row = {
            "label": display_label,
            "depth": depth,
            "n_total": int(mask.sum()),
            "n_ibd": int((sub_assign["IBD"] == 1).sum()),
            "n_control": int((sub_assign["IBD"] == 0).sum()),
        }
        for component in COMPONENTS:
            row[f"median_{component}_all_pct"] = round(float(sub_rel[component].median() * 100), 1)
            row[f"median_{component}_ibd_pct"] = round(float(ibd_rel[component].median() * 100), 1)
        rows.append(row)

    OUT.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(OUT, sep="\t", index=False)
    print(OUT)


if __name__ == "__main__":
    main()
