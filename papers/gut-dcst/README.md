# Gut DCST Paper Workspace

This workspace holds the manuscript-facing materials for the DCST gut
microbiome application paper.

## Layout

- `manuscript/`: canonical LaTeX manuscript source (`gut_application_paper.tex`)
- `assets/figures/`: manuscript-ready figure files
- `assets/tables/`: manuscript-ready table files
- `scripts/`: builders for figures, tables, and manuscript exports
- `build/`: generated PDF output and visual check artifacts
- `notes/`: action plans, handoffs, claim inventories, and figure plans
- `archive/2026-03-24-phase1/`: the phase-1 report and related archived paper files
- `archive/2026-04-03-markdown-source/`: archived Pandoc-era Markdown source

## Build Commands

From [`scripts/`](/Users/pgajer/current_projects/linf/papers/gut-dcst/scripts):

- `./build_gut_application_paper_pdf.sh`

## Source Of Record

The canonical manuscript source is
[`gut_application_paper.tex`](/Users/pgajer/current_projects/linf/papers/gut-dcst/manuscript/gut_application_paper.tex).
The older Markdown draft was archived after the workflow switched to
direct LaTeX editing for better float placement and layout control.

## Project Boundary

Large cohort outputs and analysis-scale CSV/PNG artifacts remain under
`/Users/pgajer/current_projects/gut_microbiome/outputs/`.
Only manuscript-facing figures, tables, and archived report versions should be
copied into this paper workspace.
