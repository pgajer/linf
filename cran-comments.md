## Submission

This is an update to `linf`. Version 0.3.0 consolidates hierarchy refinement in
`refine.linf.csts()`, removes three superseded public functions, and requires
the explicit dCST object structure introduced in version 0.2.0. These
intentional API changes are listed in `NEWS.md`. CRAN currently reports no
reverse dependencies for `linf`.

## Test environments

* Local macOS 26.6.1 (Apple Silicon), R-devel (2026-06-24 r90190):
  0 errors | 0 warnings | 1 note

Remote checks for version 0.3.0 have not yet been run.

## R CMD check results

0 errors | 0 warnings | 1 note

## Notes

* One NOTE records that version 0.2.0 was updated on CRAN today. Version 0.3.0
  is not being submitted immediately; this file will be refreshed before its
  eventual submission.
* The API cleanup is intentional and documented in `NEWS.md`; the package has
  no CRAN reverse dependencies.

## Additional checks

* Examples, tests, vignettes, and PDF and HTML manuals pass locally.
* The built source package contains no personal paths, debugging hooks, paper
  sources, project notes, or data-preparation directories.
