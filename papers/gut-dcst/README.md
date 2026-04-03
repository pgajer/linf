# Gut DCST Paper Workspace

This workspace holds the manuscript-facing materials for the DCST gut
microbiome application paper.

## Layout

- `manuscript/`: current paper source files (`.md`, `.tex`, `references.bib`)
- `assets/figures/`: manuscript-ready figure files
- `assets/tables/`: manuscript-ready table files
- `scripts/`: builders for figures, tables, and manuscript exports
- `build/`: generated PDF output and visual check artifacts
- `notes/`: action plans, handoffs, claim inventories, and figure plans
- `archive/2026-03-24-phase1/`: the phase-1 report and related archived paper files

## Build Commands

From [`scripts/`](/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts):

- `./build_gut_application_paper_tex.sh`
- `./build_gut_application_paper_pdf.sh`

## Project Boundary

Large cohort outputs and analysis-scale CSV/PNG artifacts remain under
`/Users/pgajer/current_projects/gut_microbiome/outputs/`.
Only manuscript-facing figures, tables, and archived report versions should be
copied into this paper workspace.
