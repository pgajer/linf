# Valencia 1k four-component hypercube embedding example

A stratified 1,000-sample subset of the Valencia 13k vaginal microbiome
training set, reduced to four selected phylotype coordinates and
L1-normalized over those coordinates. The object is bundled as a
lightweight example for the zero-aware homogeneous-coordinate hypercube
embedding in
[`linf.hypercube.embedding`](https://pgajer.github.io/linf/reference/linf.hypercube.embedding.md).

## Usage

``` r
valencia_linf_hypercube_1k
```

## Format

A list with four components:

- rel4:

  Numeric matrix (1000 x 4). Rows are samples and columns are `Li`,
  `Lc`, `Gv`, and `Bv`. Each row sums to 1 after restricting the
  original Valencia profile to the four mapped taxa.

- meta:

  Data frame (1000 rows) with columns: `sample_id` (anonymized
  bundled-data identifier), `source_row` (row number in the filtered
  Valencia 13k source), `Val_CST`, `Val_subCST`, `selected_mass`
  (original relative-abundance mass carried by the four selected taxa),
  and `dominant_component` (largest of `Li`, `Lc`, `Gv`, `Bv` after
  renormalization).

- component_map:

  Named character vector mapping `Li`, `Lc`, `Gv`, and `Bv` to the
  original Valencia taxon names: `Lactobacillus_iners`,
  `Lactobacillus_crispatus`, `Gardnerella_vaginalis`, and `BVAB1`.

- source:

  Character string documenting provenance.

## Source

Generated from the VALENCIA training data at
<https://github.com/ravel-lab/VALENCIA>. See
`data-raw/build_valencia_linf_hypercube_1k.R` and the installed
`DATA_PROVENANCE.md` file.

## Details

Rows with zero total mass in the four selected taxa are removed before
renormalization. Sampling is stratified by the dominant selected
component, using `set.seed(20261604)`. The object is not intended to
replace the full Valencia matrix; it is a compact reproducible example
for visualizing compositional projective-space coordinate charts.

## Examples

``` r
data(valencia_linf_hypercube_1k)
dim(valencia_linf_hypercube_1k$rel4)
#> [1] 1000    4
table(valencia_linf_hypercube_1k$meta$dominant_component)
#> 
#>  Li  Lc  Gv  Bv 
#> 366 321 257  56 
emb <- linf.hypercube.embedding(
  valencia_linf_hypercube_1k$rel4,
  reference = "Li"
)
dim(emb)
#> [1] 1000    3
```
