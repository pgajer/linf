# Gut DCST Paper Workspace

This workspace holds the manuscript-facing materials for the DCST gut
microbiome application paper.

## Layout

- `manuscript/`: canonical LaTeX manuscript source (`gut_application_paper.tex`)
- `assets/figures/`: manuscript-ready figure files
- `assets/tables/`: manuscript-ready table files
- `scripts/`: builders for figures, tables, and manuscript exports
- `build/`: generated PDF output and visual check artifacts
- `notes/`: action plans, handoffs, claim inventories, figure plans, and
  journal-target planning
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

The current journal-target shortlist is tracked in
[`JOURNAL_TARGETS_gut_application_paper.md`](/Users/pgajer/current_projects/linf/papers/gut-dcst/notes/JOURNAL_TARGETS_gut_application_paper.md).

The absorb-policy rerun plan is tracked in
[`ACTION_PLAN_absorb_dcst_rerun_2026-04-10.md`](/Users/pgajer/current_projects/linf/papers/gut-dcst/notes/action_plans/ACTION_PLAN_absorb_dcst_rerun_2026-04-10.md).

## Versioning Policy

- Use git and GitHub for day-to-day version history.
- Keep exactly one live manuscript source file:
  `manuscript/gut_application_paper.tex`.
- Do not create routine `_v1`, `_v2`, `_final`, or similar manuscript
  filenames.
- Reserve `archive/` for milestone snapshots only, such as:
  first full draft, preprint/submission, revision, and accepted version.
- Name milestone archive folders as `YYYY-MM-DD-<milestone>/`.
- For milestone commits, prefer lightweight tags such as
  `gut-dcst-draft-1`, `gut-dcst-submission-1`, or
  `gut-dcst-resubmission-1`.

To create a standardized archive snapshot from the current manuscript
workspace, run:

- `./archive_gut_application_paper_snapshot.sh <milestone-label>`

This creates a dated subdirectory under `archive/`, copies the current
canonical `.tex` source plus `references.bib`, copies the built PDF if it
exists, and writes a small snapshot `README.md` with the current git commit.

## Project Boundary

Large cohort outputs and analysis-scale CSV/PNG artifacts remain under
`/Users/pgajer/current_projects/gut_microbiome/outputs/`.
Only manuscript-facing figures, tables, and archived report versions should be
copied into this paper workspace.
