# Reproducing the bundled data

Run these commands from the repository root. Use R >= 4.0 (the sampling
algorithm used here is R's post-3.6 default), Python 3, and the package's
normal build dependencies. Reproduction was verified with R-devel
2026-06-24 r90190 and linf 0.3.1 on macOS arm64 on 15 September 2026.
No private input paths are needed for the VALENCIA recipe.

## VALENCIA: public CSVs to prepared matrices

The source is pinned to ravel-lab/VALENCIA revision
`8559d454387479f7155333693d854961463c3b15`. The two CSVs in
`Publication_materials/Data_and_metadata` match the preserved preparation
inputs byte for byte. [The installed manifest](../inst/DATA_MANIFEST.csv)
records URLs, SHA-256 and MD5 checksums for source and prepared files.

```sh
python3 data-raw/fetch_valencia.py build/valencia-raw
Rscript data-raw/prepare_valencia.R build/valencia-raw build/valencia-prepared
```

The preparation keeps samples with at least 3,000 reads, retains taxa observed
with at least two reads in at least 1% of those samples and at least 50 total
reads, then reapplies the 3,000-read threshold and normalizes to unit sum.
It preserves source order and aligns the supplied CST annotations. The result
is 12,881 samples by 178 taxa, with a matching annotation table and filtering
summaries. The original preserved matrices, metadata and summaries were
reproduced exactly as R objects; CSV filtering summaries were byte-identical.
Serialized RDS bytes can differ with R or compression-library versions.

## Prepared inputs to package assets

Build in a separate output directory, then compare before replacing data:

```sh
LINF_VALENCIA_ROOT=build/valencia-prepared LINF_DATA_OUTPUT=build/rebuilt-data \
  Rscript -e 'for (f in c("build_valencia2k.R", "build_valencia_linf_hypercube_1k.R", "build_valencia13k_merged_dcst_depths.R")) sys.source(file.path("data-raw", f), new.env())'
Rscript tools/audit-data.R --rebuilt=build/rebuilt-data
```

The 2,000-sample builder saves the verified sampled source-row indices in
`cst$source_row`; it does not change the established local IDs, order or values.
The four-component subset and merged assignment assets use the same source-row
numbering. Only join these keys when the prepared-input revision agrees.
The assignment builders record the generator package/R versions and upstream
revision in `params$generator`. They recompute assignments using the current
package's deterministic first-feature tie rule.

The September 2026 refresh changes two saved depth-2 assignments and ten
depth-3 assignments across ten source rows relative to the older assets.
Direct changes occur at tied candidate abundances; one downstream difference
also follows a changed parent. The input data, source-row identities, CST
annotations and depth-1 assignments are unchanged. These are refreshed
descriptive source-data results, not new independent validation.

## American Gut: shipped prepared counts to package object

The exact prepared count, taxonomy and metadata inputs are shipped in
`inst/extdata/` and pinned in the manifest. They rebuild `agp_gut` without
network access:

```sh
LINF_DATA_OUTPUT=build/rebuilt-data Rscript data-raw/build_agp_gut.R
Rscript tools/audit-data.R --rebuilt=build/rebuilt-data
```

This reproducible boundary starts at prepared counts, not raw sequencing reads.
The separate upstream PRIME/SILVA processing and original 5,000-sample fit are
not implemented by this package. `create_agp_gut_subset.py --help` documents
how to reselect the subset when those full upstream tables are available;
it is not required to reproduce the shipped package object. The sampling
rule and scientific limits are in [DATA_PROVENANCE.md](../inst/DATA_PROVENANCE.md).
Do not infer that a package-data rebuild independently validates sequence
processing or makes this stratified subset a population sample.
