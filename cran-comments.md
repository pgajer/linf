## Submission

This is a new submission of `linf`, an R package providing utilities for
L-infinity normalization, dominant-feature assignment, dominant community
state types, and landmark-profile construction for compositional data.

## Test environments

* Local macOS 26.3.1 (Apple Silicon), R 4.5.2:
  0 errors | 0 warnings | 1 note
* Local macOS 26.3.1 (Apple Silicon), R-devel (2026-06-24 r90190):
  0 errors | 0 warnings | 1 note
* Win-builder, Windows Server 2022, R 4.6.1:
  0 errors | 0 warnings | 1 note
* Win-builder, Windows Server 2022, R-devel (2026-07-25 r90301):
  0 errors | 0 warnings | 1 note
* Win-builder, Windows Server 2022, R 4.5.3:
  0 errors | 0 warnings | 1 note
* R-hub, Linux, Windows, and macOS, R-devel:
  0 errors | 0 warnings | 0 notes on each platform
* GitHub Actions:
  Ubuntu (R-oldrel, R-release, and R-devel), Windows (R-release), and
  macOS (R-release): 0 errors | 0 warnings | 0 notes on each job

## R CMD check results

0 errors | 0 warnings | 1 note

## Notes

* `New submission`: this is the first CRAN submission of `linf` and the only
  NOTE in both local checks.
* The isolated local R-devel library did not provide the optional suggested
  package `plotly`; this was reported as INFO, not as a NOTE.

## Additional checks

* `urlchecker::url_check()` reports that all URLs are correct.
* Examples, tests, vignettes, and PDF and HTML manuals pass locally and on
  Win-builder.
* The bundled AGP-derived subset is reproducible from an explicit,
  deterministic, seed-controlled selection script.
