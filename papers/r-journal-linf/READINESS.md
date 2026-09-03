# R Journal readiness

## Package correction gate (2026-09-03)

The archived article remains reproducible with its pinned CRAN 0.3.0 source,
but is not ready for journal submission. The 0.3.1 transfer correction removes
one out-of-hierarchy fallback assignment in the held-out example: depth-1
coverage becomes 499/500, while depth-2 coverage remains 486/500. Depth-1
concordance among assigned samples remains 93.8% after rounding. After 0.3.1
is public on CRAN, update the source pin, citation, coverage text and assertions,
then regenerate and audit the article and archive. The successful 0.3.0 build
below is reproduction evidence, not clearance of this correctness issue.

## Verified state (2026-09-02)

- The article builds from the checksum-pinned CRAN `linf` 0.3.0 source and
  produces a visually reviewed eight-page PDF plus a self-contained HTML
  article.
- All applicable `rjtools` checks pass, including file structure, title and
  section case, abstract format, spelling, package availability, bibliography,
  and the fixed submission date.
- The citation audit maps all 10 cited keys to verified source links and
  supporting claims. The `linfPackage` row was reverified against the published
  CRAN 0.3.0 source and manual, and its version is checked against the source
  pin and BibTeX entry.
- `R CMD check --as-cran` under the locally installed R-devel completes with
  no errors, warnings, or notes with CRAN incoming-upload checks disabled for
  this already-published version. Examples, tests, vignettes, and HTML/PDF
  manuals are included. An unmodified incoming check warns only that 0.3.0 is
  already on CRAN; this is not a journal-submission defect.
- The review archive includes the Makefile, audit scripts, exact 0.3.0 package
  source tarball, manuscript sources, motivating letter, LaTeX figure files,
  rendered outputs, and SHA-256 manifests.
  Its automated test clears user R-library/profile variables, installs the
  bundled tarball into a fresh library, removes extracted generated outputs,
  and reruns the article audit in a path containing spaces. Missing build
  dependencies are installed before the timed reproduction test.
- The benchmark table's percent sign is escaped for LaTeX, and tables are
  placed after their introductions. The readiness scan rejects escaped table
  markup, missing equation/table labels, and missing HTML accessibility labels.
- Transfer coverage is reported separately at depths 1 and 2, and absorb-policy
  labels are not interpreted as necessarily naming a sample's largest feature.
  Benchmark storage is explicitly input-object size, not peak workflow memory.
- The compositions 2.0-9 citation year is corrected to 2025 against the official
  CRAN manual; its source-verified year is checked against the bibliography and
  visible verification record. Figure 2 has headroom for all count labels.
- Every figure referenced by the generated TeX must exist and be included in
  the archive; a final PDF render preserves the supporting figure files.
- The printed transfer example includes the split and reference fitting;
  both figure captions explain their annotations. The archive records the
  PDF benchmark's unrounded results, session information, package versions,
  and basic machine information. The Matrix citation matches version 1.7-5
  used in this draft; other build environments may use different versions.
- The pkgdown website is live at <https://pgajer.github.io/linf/> (HTTP 200).

## Automated gates

- `make dependencies`: installs missing article build dependencies.
- `make render`: installs the pinned CRAN `linf` source into an isolated library
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

## Public-version gate

The article describes `linf` 0.3.0, published on CRAN on 2026-08-31. The
source pin, isolated installation, and public CRAN version agree for that
historical build. The package correction gate above now requires an updated
public source version and article rebuild before submission.

## Human gates

The author must approve authorship metadata, the AI-use disclosure, all
publication-facing prose and results, the motivating letter, and the exact
submission archive; and must confirm that the article is not published,
submitted, or under review elsewhere. The article date must also be changed to
the actual submission date if submission occurs after the current draft date,
2026-09-02. These human approvals have not been inferred from the package's
CRAN acceptance.
