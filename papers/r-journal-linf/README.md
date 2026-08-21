# `linf` R Journal paper

This directory is the canonical workspace for the software-focused R Journal
article describing the `linf` package. The mathematical development remains in
the separate `Linf_paper` arXiv workspace; this article cites that work and
focuses on package design, workflows, implementation, reproducibility, and
limitations.

The canonical manuscript source is `linf.Rmd`. Generated article wrappers,
HTML, PDF, figures, and submission archives are build artifacts.

## Build

From this directory:

```sh
make render
make citation-check
make rj-check
make audit
make draft-bundle
```

`make render` installs the current package source into an isolated library and
then evaluates the article. `make submission-audit` additionally requires the
manuscript's target package version to be public on CRAN and runs the package's
full clean CRAN-style check.

Final rendered artifacts are written to `output/pdf/` and `output/html/`.

The archive produced by `make draft-bundle` contains this Makefile, the audit
scripts, the rendered article and motivating letter, and an installable source
tarball for the exact `linf` version used by the article. After extracting the
archive, a reviewer can run `make audit`; the bundled tarball is installed into
the archive's local `library/` before every result is regenerated. The
`version-check`, `package-check`, and `submission-audit` targets are intended
for the source repository because they additionally inspect CRAN and the full
package checkout.
