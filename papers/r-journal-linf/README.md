# `linf` R Journal paper

This directory is the canonical workspace for the software-focused R Journal
article describing the `linf` package. The mathematical development remains in
the separate `Linf_paper` arXiv workspace; this article cites that work and
focuses on package design, workflows, implementation, reproducibility, and
limitations.

The canonical manuscript source is `linf.Rmd`. Generated article wrappers,
HTML, PDF, figures, and submission archives are build artifacts.

## Build

The build requires R, Python 3.9 or newer, Make, Pandoc, a LaTeX installation,
and Poppler (`pdfinfo`, `pdffonts`, and `pdftotext`). Install the R packages
listed in `_Rpackages.txt`; the build itself supplies the pinned `linf`
version. Internet access is needed for the first source download and for the
journal's package-availability checks. A current HTML Tidy is recommended for
the optional CRAN-style package check.

From this directory, or from an extracted review archive:

```sh
make dependencies
make render
make citation-check
make rj-check
make audit
make draft-bundle
```

`make render` obtains the exact CRAN source specified in `package-source.json`,
verifies its SHA-256 checksum, installs it into a local library, and evaluates
the article. An existing bundled tarball is verified and reused without a
download. No package checkout or maintainer library is required.
`make dependencies` installs any missing article dependencies into the local
library; this one-time setup is separate from article reproduction timing.
`make submission-audit` additionally checks public CRAN availability and runs
`R CMD check --as-cran` on the same pinned source tarball.
Only CRAN's incoming-upload checks are disabled for this published-source
check: they would reject uploading version 0.3.0 again because it is already
on CRAN. Examples, tests, vignettes, and manual checks remain enabled, and the
target fails unless the check log reports `Status: OK`.

Final rendered artifacts are written to `output/pdf/` and `output/html/`.

The archive produced by `make draft-bundle` includes this Makefile, all build
and audit scripts, the article sources and rendered outputs, the LaTeX wrapper
and figure files, the motivating letter, and the pinned CRAN source tarball. It also
includes `SHA256SUMS` and a source/output manifest. Bundle creation rejects
source/output drift and missing files referenced by LaTeX `includegraphics`.
An explicit final PDF render preserves these figure files after the combined
HTML/PDF render's cleanup. Bundle creation tests `make audit` after removing the extracted
generated outputs and clearing user R-library/profile settings. Reviewers
can run the same commands in the extracted archive. Numerical results are
recomputed; benchmark timings can vary across machines and runs.

The archive's `build/benchmark-results.csv` contains the unrounded PDF benchmark
results. `build/session-info.txt` records `sessionInfo()`, including BLAS/LAPACK,
and `build/benchmark-environment.json` records loaded package versions, R,
operating system, architecture, logical cores, and CPU/memory details when
available. Hostnames and user names are not collected. These records describe
the actual build; dependencies other than linf are not version-locked.
