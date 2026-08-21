# R Journal readiness

## Verified state (2026-08-21)

- The article builds from the local `linf` 0.3.0 source and
  produces a visually reviewed eight-page PDF plus a self-contained HTML
  article.
- All applicable `rjtools` checks pass, including file structure, title and
  section case, abstract format, spelling, package availability, bibliography,
  and the fixed submission date.
- The citation audit maps all 10 cited keys to verified source links and
  supporting claims, and identifies 0.3.0 as the manuscript's development
  version while retaining a separate public-CRAN gate.
- `R CMD check --as-cran` under the locally installed R-devel completes with
  no errors or warnings and one expected timing NOTE, including examples,
  tests, vignettes, and HTML/PDF manuals.
- The review archive includes the Makefile, audit scripts, exact 0.3.0 package
  source tarball, manuscript sources, motivating letter, and rendered outputs.
  Its automated test clears user R-library/profile variables, installs the
  bundled tarball into a fresh library, and reruns the article audit after
  extraction.

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

The article describes `linf` 0.3.0. Submission is blocked until that version is
publicly available on CRAN. Drafting and technical review can proceed against
the isolated installation of the local 0.3.0 source. The automated gate
currently finds CRAN version 0.2.0 and therefore fails as intended.

## Human gates

The author must approve authorship metadata, the AI-use disclosure, all
publication-facing prose and results, the motivating letter, and the exact
submission archive; and must confirm that the article is not published,
submitted, or under review elsewhere. The article date must also be changed to
the actual submission date if submission occurs after 2026-08-21.
