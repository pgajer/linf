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
* All 230 test assertions pass, with no failures, warnings, or skips.
* Examples, both vignettes (including rebuilding), and PDF/HTML manuals pass.
* GitHub Actions: macOS, Windows and Ubuntu R 4.6.1 release, and Ubuntu
  R 4.5.3 oldrelease: all Status: OK; 230 passing assertions on each platform.
* R-hub Windows, R-devel 2026-09-02 r90474 ucrt: Status: OK;
  230 passing assertions.

GitHub Ubuntu devel and R-hub Linux/macOS are still running. The exact source
tarball has been sent to Win-builder release/devel/oldrelease; result logs are
pending. This is a preparation record, not a statement that every requested
pre-submission check is complete.
