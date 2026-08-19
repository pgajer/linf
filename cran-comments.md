## Submission

This is an update to `linf`. Version 0.2.0 adopts the dominance-sample-set and
dominance-lineage terminology used in the methods paper. This intentionally
changes the public function and result-field names; the changes are listed in
`NEWS.md`. CRAN currently reports no reverse dependencies for `linf`.

## Test environments

* Local macOS 26.6.1 (Apple Silicon), R-devel (2026-06-24 r90190):
  0 errors | 0 warnings | 1 note
* Win-builder, Windows Server 2022, R 4.6.1:
  0 errors | 0 warnings | 0 notes
* Win-builder, Windows Server 2022, R-devel (2026-08-17 r90424):
  0 errors | 0 warnings | 0 notes
* Win-builder, Windows Server 2022, R 4.5.3:
  0 errors | 0 warnings | 0 notes
* R-hub, Linux, Windows, and macOS, R-devel:
  0 errors | 0 warnings | 0 notes on each platform
* GitHub Actions, Ubuntu release/devel/oldrel, Windows release, and macOS
  release: 0 errors | 0 warnings | 0 notes on each job

## R CMD check results

0 errors | 0 warnings | 1 note

## Notes

* The local NOTE is environmental: HTML-manual validation was skipped because
  the installed HTML Tidy is not recent enough. Win-builder reports no NOTE.
* The compatibility-breaking terminology changes are intentional and are
  documented in `NEWS.md`; the package has no CRAN reverse dependencies.

## Additional checks

* Examples, tests, vignettes, and PDF and HTML manuals pass locally.
* The built source package contains no personal paths, debugging hooks, paper
  sources, project notes, or data-preparation directories.
