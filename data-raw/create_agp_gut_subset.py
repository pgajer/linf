#!/usr/bin/env python3
"""Create the deterministic 766-sample AGP demonstration subset.

Selection is based only on retained 5,000-sample DCST assignments:

1. Keep every sample in four uncommon demonstration DCSTs.
2. Exclude Eukaryota and Unassigned depth-1 labels.
3. Fill the remaining slots by a seeded simple random sample from the
   eligible background.

Phenotype fields are joined only after sample membership has been fixed and
therefore cannot influence selection.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import io
import math
import os
import random
from collections import Counter
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
TARGET_DCSTS = {
    "Akkermansia sp.",
    "Pasteurellaceae",
    "Prevotella_7 sp.",
    "Staphylococcus sp.",
}
EXCLUDED_DCSTS = {"Eukaryota", "Unassigned"}
PHENOTYPE_PATTERNS = {
    "IBS": "Irritable bowel syndrome (IBS)",
    "IBD": "Inflammatory bowel disease (IBD)",
    "Diabetes": "Diabetes",
    "Autoimmune": "Autoimmune",
    "Seasonal_allergies": "Seasonal allergies",
    "Migraine": "Migraine",
    "Acid_reflux": "Acid reflux",
    "Lung_disease": "Lung disease",
    "Cardiovascular_disease": "Cardiovascular disease",
    "Skin_condition": "Skin condition",
}
META_FIELDS = [
    "Run",
    "dcst_depth1",
    "dcst_depth2",
    *PHENOTYPE_PATTERNS,
    "Obesity",
    "BMI",
]


def required_path(value: str | None, option: str) -> Path:
    if not value:
        raise SystemExit(
            f"{option} is required (pass it explicitly or set its documented "
            "environment variable)."
        )
    path = Path(value).expanduser().resolve()
    if not path.is_file():
        raise SystemExit(f"{option} does not exist: {path}")
    return path


def read_assignments(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        rows = list(reader)
    required = {"Run", "dcst_depth1", "dcst_depth2"}
    missing = required.difference(reader.fieldnames or [])
    if missing:
        raise SystemExit(f"DCST assignments lack columns: {sorted(missing)}")
    run_ids = [row["Run"] for row in rows]
    if len(run_ids) != len(set(run_ids)):
        raise SystemExit("DCST assignments contain duplicated Run values.")
    return rows


def select_runs(
    assignments: list[dict[str, str]], n_samples: int, seed: int
) -> tuple[list[str], dict[str, str]]:
    target_rows = [
        row for row in assignments if row["dcst_depth1"] in TARGET_DCSTS
    ]
    background_rows = [
        row
        for row in assignments
        if row["dcst_depth1"] not in TARGET_DCSTS | EXCLUDED_DCSTS
    ]
    if len(target_rows) > n_samples:
        raise SystemExit(
            f"Target DCSTs contain {len(target_rows)} rows, exceeding "
            f"n_samples={n_samples}."
        )
    n_background = n_samples - len(target_rows)
    if n_background > len(background_rows):
        raise SystemExit(
            f"Only {len(background_rows)} eligible background rows are "
            f"available for {n_background} slots."
        )

    rng = random.Random(seed)
    sampled_background = set(
        rng.sample(sorted(row["Run"] for row in background_rows), n_background)
    )
    selected = {
        row["Run"]: (
            "all_target_dcst"
            if row["dcst_depth1"] in TARGET_DCSTS
            else "seeded_random_background"
        )
        for row in assignments
        if row["dcst_depth1"] in TARGET_DCSTS
        or row["Run"] in sampled_background
    }
    ordered_runs = [row["Run"] for row in assignments if row["Run"] in selected]
    if len(ordered_runs) != n_samples:
        raise SystemExit(
            f"Internal selection error: selected {len(ordered_runs)} rows."
        )
    return ordered_runs, selected


def read_agp_metadata(path: Path) -> dict[str, dict[str, str]]:
    with gzip.open(path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {"Run", "BioProject", "Phenotype", "Host_BMI"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise SystemExit(f"Metadata lacks columns: {sorted(missing)}")
        rows = {
            row["Run"]: row
            for row in reader
            if row.get("BioProject") == "PRJEB11419"
        }
    return rows


def phenotype_row(
    assignment: dict[str, str], source: dict[str, str]
) -> dict[str, str | int]:
    phenotype = (source.get("Phenotype") or "").lower()
    output: dict[str, str | int] = {
        "Run": assignment["Run"],
        "dcst_depth1": assignment["dcst_depth1"],
        "dcst_depth2": assignment["dcst_depth2"],
    }
    for name, pattern in PHENOTYPE_PATTERNS.items():
        output[name] = int(pattern.lower() in phenotype)

    bmi_text = (source.get("Host_BMI") or "").strip()
    try:
        bmi = float(bmi_text)
    except ValueError:
        bmi = math.nan
    output["Obesity"] = int(not math.isnan(bmi) and bmi >= 30)
    output["BMI"] = "" if math.isnan(bmi) else bmi_text
    return output


def extract_counts(
    abundance_path: Path, run_ids: list[str], taxa: list[str]
) -> dict[str, list[str]]:
    wanted_runs = set(run_ids)
    selected: dict[str, list[str]] = {}
    with gzip.open(abundance_path, "rt", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        header = next(reader)
        if not header or header[0] != "Run":
            raise SystemExit("Abundance table must begin with a Run column.")
        positions = {name: index for index, name in enumerate(header)}
        missing_taxa = [name for name in taxa if name not in positions]
        if missing_taxa:
            raise SystemExit(
                f"{len(missing_taxa)} bundled taxa are absent from abundance table."
            )
        taxa_indices = [positions[name] for name in taxa]
        for row in reader:
            run_id = row[0]
            if run_id in wanted_runs:
                selected[run_id] = [row[index] for index in taxa_indices]
    missing_runs = [run_id for run_id in run_ids if run_id not in selected]
    if missing_runs:
        raise SystemExit(
            f"{len(missing_runs)} selected runs are absent from abundance table."
        )
    return selected


def write_gzip_csv(path: Path, rows: list[list[str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0) as zipped:
            with io.TextIOWrapper(zipped, encoding="utf-8", newline="") as text:
                writer = csv.writer(text, lineterminator="\n")
                writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--abundance",
        default=os.getenv("LINF_AGP_ABUNDANCE"),
        help="Full PRIME SILVA abundance CSV.gz (LINF_AGP_ABUNDANCE).",
    )
    parser.add_argument(
        "--metadata",
        default=os.getenv("LINF_AGP_METADATA"),
        help="PRIME sample metadata CSV.gz (LINF_AGP_METADATA).",
    )
    parser.add_argument(
        "--dcst-assignments",
        default=os.getenv("LINF_AGP_DCST_ASSIGNMENTS"),
        help=(
            "Retained 5,000-sample dcst_assignments.csv "
            "(LINF_AGP_DCST_ASSIGNMENTS)."
        ),
    )
    parser.add_argument(
        "--taxa",
        default=str(REPO_ROOT / "inst" / "extdata" / "agp_gut_taxa.txt"),
        help="Ordered taxon list retained from the 5,000-sample analysis.",
    )
    parser.add_argument(
        "--output-counts",
        default=str(
            REPO_ROOT / "inst" / "extdata" / "agp_gut_counts.csv.gz"
        ),
        help="Output count matrix.",
    )
    parser.add_argument(
        "--output-metadata",
        default=str(REPO_ROOT / "inst" / "extdata" / "agp_gut_meta.csv"),
        help="Output metadata and DCST manifest.",
    )
    parser.add_argument("--n-samples", type=int, default=766)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    abundance = required_path(args.abundance, "--abundance")
    metadata_path = required_path(args.metadata, "--metadata")
    assignments_path = required_path(
        args.dcst_assignments, "--dcst-assignments"
    )
    taxa_path = required_path(args.taxa, "--taxa")

    assignments = read_assignments(assignments_path)
    assignment_by_run = {row["Run"]: row for row in assignments}
    run_ids, selection_reason = select_runs(
        assignments, args.n_samples, args.seed
    )
    metadata = read_agp_metadata(metadata_path)
    missing_metadata = [run_id for run_id in run_ids if run_id not in metadata]
    if missing_metadata:
        raise SystemExit(
            f"{len(missing_metadata)} selected runs lack AGP metadata."
        )
    taxa = [
        line
        for line in taxa_path.read_text(encoding="utf-8").splitlines()
        if line
    ]
    counts = extract_counts(abundance, run_ids, taxa)

    count_rows = [["Run", *taxa]]
    count_rows.extend([[run_id, *counts[run_id]] for run_id in run_ids])
    write_gzip_csv(Path(args.output_counts).expanduser().resolve(), count_rows)

    output_metadata = Path(args.output_metadata).expanduser().resolve()
    output_metadata.parent.mkdir(parents=True, exist_ok=True)
    with output_metadata.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=[*META_FIELDS, "selection_reason"], lineterminator="\n"
        )
        writer.writeheader()
        for run_id in run_ids:
            row = phenotype_row(assignment_by_run[run_id], metadata[run_id])
            row["selection_reason"] = selection_reason[run_id]
            writer.writerow(row)

    label_counts = Counter(
        assignment_by_run[run_id]["dcst_depth1"] for run_id in run_ids
    )
    reason_counts = Counter(selection_reason.values())
    print(
        f"Selected {len(run_ids)} rows from {len(assignments)} assignments "
        f"with seed {args.seed}."
    )
    print(f"Selection reasons: {dict(sorted(reason_counts.items()))}")
    for label, count in label_counts.most_common():
        print(f"{label}: {count}")
    print(f"Wrote counts: {Path(args.output_counts).expanduser().resolve()}")
    print(f"Wrote metadata: {output_metadata}")


if __name__ == "__main__":
    main()
