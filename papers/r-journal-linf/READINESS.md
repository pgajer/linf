# R Journal readiness

## Verified state (2026-08-20)

- The article builds from the local `linf` 0.2.0 source in 26 seconds and
  produces a visually reviewed seven-page PDF plus a self-contained HTML
  article.
- All applicable `rjtools` checks pass, including file structure, title and
  section case, abstract format, spelling, package availability, bibliography,
  and the fixed submission date.
- The citation audit maps all 10 cited keys to verified source links and
  supporting claims.
- `R CMD check --as-cran` under the locally installed R-devel completes with
  `Status: OK`, including examples, tests, vignettes, and HTML/PDF manuals.
- The minimal 12-file review archive rebuilds successfully after extraction.

## Automated gates

- `make render`: installs the current `linf` source into an isolated library
  and evaluates the article in the official `rjtools` format.
- `make citation-check`: requires one verified evidence row for every cited
  key and rejects unresolved bibliography entries.
- `make rj-check`: runs the applicable `rjtools` initial checks.
- `make audit`: scans the public source and rendered artifacts for private
  paths, prompt traces, unresolved markers, undefined references, malformed
  PDFs, page-count excess, and unembedded fonts.
- `make draft-bundle`: creates and clean-builds a minimal review archive.
- `make submission-audit`: adds the public CRAN version gate and a full clean
  CRAN-style package check.

## Current submission blocker

The article describes `linf` 0.2.0. Submission is blocked until that version is
publicly available on CRAN. Drafting and technical review can proceed against
the isolated installation of the local 0.2.0 source. The automated gate
currently finds CRAN version 0.1.0 and therefore fails as intended.

## Human gates

The author must approve authorship metadata, the AI-use disclosure, all
publication-facing prose and results, the motivating letter, and the exact
submission archive. The article date must also be changed to the actual
submission date if submission occurs after 2026-08-20.
