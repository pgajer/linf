## Update

Version 0.3.1 fixes frozen transfer with distinct feature IDs and display labels,
prevents transfer outside the fitted hierarchy, and makes dense and sparse
absorption agree on ties. It preserves taxonomy identifiers containing the
lineage separator during refinement, transfer and landmark lookup, validates feature-metadata
lengths, and clarifies zero-row and tolerance behavior.

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
* All 476 test assertions pass, with no failures, warnings or skips.
* Examples, all four installed vignettes and their rebuilding, and the PDF
  and HTML manuals pass. The installed index includes HTML, Rmd and R sources.
* Regression tests cover separator-containing IDs, suffix collisions, terminal
  lineages, custom rare labels, dense/sparse matrices and all policy views.

The earlier September 3 platform checks apply to an older candidate. Fresh Windows,
Linux, macOS and Win-builder checks are being run for this candidate; results
will be recorded here before release coordination.

## Preparation status

Not submitted. CRAN currently publishes 0.3.0. No linf entry was found in the
public pending/inspect queues or a 0.3.1 submission in the connected mail
search on 2026-09-14.

The landmark lookup defect documented during the vignette review is fixed.
Feature identities are recovered from complete fitted hierarchy transitions,
so separators inside IDs do not cause missing targets or suffix collisions.
