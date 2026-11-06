## Update

This is linf 0.4.0, an update to CRAN version 0.3.0.

The release adds `linf.udcsts()` to group samples by their leading features
regardless of order, and `linf.udcst.graph()` to connect two-feature groups
supported by observed three-feature groups. It fixes feature matching and
assignment of new samples to existing groups, makes choices between equal
values consistent for ordinary and sparse matrices, and improves numerical
stability and input validation. Older grouping results are updated on use when
their groups can be identified unambiguously; otherwise they must be recreated.

Two installed guides, expanded tutorials, data provenance and corrected saved
tie assignments accompany the changes. Abundance matrices and sample identities
are unchanged. The license, maintainer and dependencies are unchanged.

## Checks (2026-11-06)

* macOS Tahoe 26.6.1 arm64, R-devel 2026-06-24 r90190: full incoming-enabled
  `R CMD check --as-cran`, with `R_MAKEVARS_USER=/dev/null` and HTML Tidy 5.8.0;
  0 errors, 0 warnings, 0 notes. All 873 test assertions pass with no failures,
  warnings or skips. Examples, all four vignettes and their rebuilding, and
  PDF and HTML manuals pass.
* Win-builder, Windows Server 2022 x64, R-devel 2026-11-05 r90840 ucrt:
  0 errors, 0 warnings, 1 note. All 873 test assertions pass with no failures,
  warnings or skips; examples, vignette rebuilding and both manuals pass.
  The service checked the same source archive as the local check.
* Bundled-data checksums, source/installed guide checks and the pkgdown build
  pass. Installed-guide links to unpublished website articles now point to
  available, explicitly labelled article sources.

## Note

Win-builder flags "Gajer" as possibly misspelled in DESCRIPTION. This is the
correctly spelled surname in the method citation.

## Downstream checks

CRAN lists no reverse dependencies across Depends, Imports, LinkingTo,
Suggests and Enhances (checked 2026-11-06).

The local optional consumer gflowui has existing dependency-declaration and
test failures. Earlier full checks reproduced the same failures with the
reviewed 0.3.1 development baseline and the preceding 0.4.0 archive. The final
archive changes only documentation and the build timestamp; package code,
data and tests are unchanged. The focused linf integration suite was rerun
against the final archive and passes all 54 assertions.
