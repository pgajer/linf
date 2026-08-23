# Repository Instructions

## Scope

- This repository is the source of truth for the `linf` R package, its
  vignettes, bundled demo data, and paper-facing gut DCST materials that
  live in-repo.

## Preferred Skills

- Prefer `$r-package-qa` for package QA, CRAN-style checks,
  documentation drift, and release-readiness work.
- Prefer `$manuscript-build-review` for paper and manuscript work under
  `papers/`.

## Package Source

- Core package code lives in `R/`, `man/`, `tests/testthat/`,
  `vignettes/`, and `inst/`.
- Demo or bundled data sources live under `data/` and `data-raw/`.
- Notes and manuscript-side material live under `notes/` and `papers/`.

## R Package Hygiene

- Edit roxygen comments in source files, not generated `.Rd` files.
- Edit `README.Rmd`, not generated `README.md` or `README.html`.
- Prefer Makefile targets over ad hoc package commands:
  - `make document`
  - `make check-fast`
  - `make check`

## Generated Files

- Do not hand-edit generated files such as `NAMESPACE`, `man/*.Rd`,
  generated `README.md`, tarballs, or `.Rcheck/` outputs unless the task
  explicitly targets generated artifacts.
- Keep local build artifacts and logs out of functional commits unless
  the task is release-specific.

## Data and Writing Safety

- Treat `data/` as curated package data and `data-raw/` as
  source-preparation material; do not rewrite them casually.
- Preserve documented public API names, vignette intent, and manuscript
  claims unless the user explicitly requests a change.

## Private Agent Work Products

- Store internal agent-only work products under
  `~/.codex/private/linf/`, not in the repository. This includes
  internal audits, agent-to-agent handoffs, intermediate rewrites,
  working copies of reviewer reports used for agent tasks, historical
  prompts, and generated review diffs that are not intended as package,
  manuscript, reproducibility, or submission artifacts.
- Organize private material first by workstream and then, when useful,
  by artifact type. Use clearly named workstream directories with
  subdirectories such as `audits/`, `handoffs/`, `drafts/`, `prompts/`,
  and `diffs/`.
- Maintain a `README.md` in each workstream directory identifying every
  file’s origin, former repository location, purpose, and possible
  future disposition.
- Keep formal and publication-facing assets in the repository. Do not
  move source code, tests, package documentation, manuscript source,
  bibliography, figures, rendering tools, citation-verification
  evidence, reproducibility inputs or scripts, checksums, provenance
  records, or formal submission files into the private tree.
- Treat draft responses, internal referee simulations, and agent working
  copies of received reports as private. If a response-to-reviewers
  document becomes part of an actual submission, copy its finalized
  version into the appropriate repository submission bundle.
- Do not make repository builds, tests, manuscript renders, or
  validation workflows depend on files under the private tree.
- When retiring a tracked internal file from the repository, preserve
  its existing Git history through the normal repository deletion and
  record its private destination in the workstream README.
- The private directory is not a credentials store. Never place
  passwords, access tokens, private keys, or other authentication
  secrets there.
