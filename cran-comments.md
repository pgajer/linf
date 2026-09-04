## Update

Version 0.3.1 fixes frozen transfer with distinct feature IDs and display labels,
prevents transfer outside the fitted hierarchy, and makes dense and sparse
absorption agree on ties. It also fixes taxonomy identifiers containing the
lineage separator, validates feature-metadata lengths, and clarifies zero-row
and tolerance behavior. Regression tests cover these corrections.

The license, maintainer, exports, dependencies, and bundled data are unchanged.
CRAN lists no reverse dependencies (checked 2026-09-03).

## Checks (2026-09-03)

* macOS 26.6.1 arm64, R-devel 2026-06-24 r90190, full incoming-enabled
  `R CMD check --as-cran` with clean user build flags and current HTML Tidy:
  0 errors, 0 warnings, 1 NOTE: "Days since last update: 3".
* All 234 test assertions pass, with no failures, warnings, or skips.
* Examples, both vignettes (including rebuilding), and PDF/HTML manuals pass.
* GitHub Actions: macOS and Windows release, Ubuntu release/devel/oldrel-1:
  all five checks report Status: OK.
* R-hub: Linux, Windows, and macOS R-devel: all three checks report Status: OK.
  R-devel revisions are r90473, r90474, and r90478, respectively.
* Win-builder release (R 4.6.1 ucrt) and oldrelease (R 4.5.3 ucrt):
  each has 0 errors, 0 warnings, 1 NOTE.

The Win-builder NOTE reports "Days since last update: 4" and flags "Gajer"
as a possible misspelling. Gajer is the correctly spelled author surname.
The update-interval difference from the local check reflects the check dates.
The CI checks use `--no-manual --as-cran`; the local full check includes manuals
and incoming feasibility checks. All completed platform suites pass 234 assertions.

The final Win-builder devel result is pending; this preparation record
will be finalized after its complete logs are reviewed.
