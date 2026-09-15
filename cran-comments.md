## Update

Version 0.3.1 fixes frozen transfer with distinct feature IDs and display labels,
prevents transfer outside the fitted hierarchy, and makes dense and sparse
absorption agree on ties. It preserves taxonomy identifiers containing the
lineage separator during refinement, transfer and landmark lookup, validates feature-metadata
lengths, and clarifies zero-row and tolerance behavior.

This candidate also separates fitted node identity from displayed paths, including
literal-path and rare-sentinel collisions; validates count filtering and empty
results; and reports assigned, unassigned and rare sample counts explicitly.
Unambiguous legacy fits upgrade on use; ambiguous legacy fits require rebuilding.

This candidate also adds two installed vignettes: a task-oriented function guide
and an example-dataset guide. Package help, README and the existing tutorials
now distinguish fitting from frozen transfer, observed landmarks from averaged
profiles, and descriptive source-data comparisons from independent validation.

The license, maintainer, function exports, dependencies and bundled data are
unchanged. CRAN lists no reverse dependencies across Depends, Imports, LinkingTo,
Suggests and Enhances (checked 2026-09-14).

## Checks (2026-09-15)

* macOS Tahoe 26.6.1 arm64; R-devel 2026-06-24 r90190.
* Full incoming-enabled `R CMD check --as-cran` of the current 0.3.1 tarball,
  with `R_MAKEVARS_USER=/dev/null` and HTML Tidy 5.8.0:
  0 errors, 0 warnings, 0 NOTEs.
* All 628 test assertions pass, with no failures, warnings or skips.
* Examples, all four installed vignettes and their rebuilding, and the PDF
  and HTML manuals pass. The installed index includes HTML, Rmd and R sources.
* The installed package overview, both help aliases, all four vignette entries
  and local inter-vignette links pass the installed guide audit.
* Regression tests cover separator-containing IDs, suffix collisions, terminal
  lineages, custom rare labels, dense/sparse matrices and all policy views.

## Current platform checks

Checked implementation commit: f537826 (2026-09-15).

* GitHub Actions: R 4.6.1 on Ubuntu, Windows and macOS; R 4.5.3 on Ubuntu;
  R-devel 2026-09-14 r90539 on Ubuntu. All five report Status: OK, with
  0 errors, 0 warnings, 0 NOTEs and 628 passing assertions.
* These jobs use `--as-cran --no-manual`, with incoming checks disabled.
  The full local check above additionally covers incoming checks and manuals.
* The pkgdown workflow builds the website successfully. Its deployment step
  is skipped on push; no website publication was performed.

Workflow evidence:
https://github.com/pgajer/linf/actions/runs/34992320125
https://github.com/pgajer/linf/actions/runs/34992320081

The following R-hub and Win-builder results are historical; those services
have not checked the new node representation.

## Earlier platform checks (2026-09-15)

These checks used the earlier package sources in commit 12422e8.

* GitHub Actions: R 4.6.1 on Ubuntu, Windows and macOS; R 4.5.3 on Ubuntu;
  R-devel 2026-09-14 r90539 on Ubuntu. All five report 0 errors, 0 warnings,
  0 NOTEs and 476 passing assertions.
* R-hub: R-devel 2026-09-14 r90539 on Linux, Windows and macOS. All three
  report 0 errors, 0 warnings, 0 NOTEs and 476 passing assertions. The first
  macOS attempt was manually interrupted during dependency compilation; the
  completed retry used unchanged package sources and workflow configuration.
* The GitHub and R-hub jobs use `--as-cran --no-manual`, with incoming checks
  disabled. The full local and Win-builder checks additionally cover incoming
  checks and manuals.
* Win-builder: R 4.6.1, R 4.5.3 and R-devel 2026-09-14 r90539, all on Windows
  Server 2022. Each reports 0 errors, 0 warnings and 1 NOTE, with all 476
  assertions passing. The sole NOTE flags "Gajer" as possibly misspelled in
  DESCRIPTION; this is the correctly spelled surname in the method citation.
  All three received the same locally checked source archive.

Workflow evidence:
https://github.com/pgajer/linf/actions/runs/34982143909
https://github.com/pgajer/linf/actions/runs/34982081022

## Preparation status

Not submitted. CRAN currently publishes 0.3.0. No linf entry was found in the
public pending/inspect queues or a 0.3.1 submission in the connected mail
search on 2026-09-14.

The landmark lookup defect documented during the vignette review is fixed.
Feature identities are read from explicit fitted node paths, so separators
inside IDs do not cause missing targets, suffix collisions or merged paths.
Website preparation is complete; publication remains separately coordinated.
