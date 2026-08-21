# Valencia 13k merged depth-3 dCST assignments

A lightweight bundled assignment asset containing merged depth-3
Dominant community state type (dCST) labels for the full Valencia 13k
vaginal microbiome training set.

## Usage

``` r
valencia13k_dcst_depth3_merged
```

## Format

A list with five components:

- assignments:

  Data frame with one row per source sample and columns `sample_id`
  (anonymized bundled-data identifier), `source_row` (row number in the
  filtered Valencia 13k source), `Val_CST`, `Val_subCST`, `dcst_depth1`,
  `dcst_depth2`, and `dcst_depth3`.

- summaries:

  Named list of per-depth summary tables. Each table contains `depth`,
  `dcst_label`, `n`, `prop`, and `path_length`.

- feature_labels:

  Character vector of source taxon labels.

- params:

  List recording the construction parameters.

- source:

  Character string documenting provenance.

## Source

Generated from the VALENCIA training data at
<https://github.com/ravel-lab/VALENCIA>. See
`data-raw/build_valencia13k_merged_dcst_depths.R` and the installed
`DATA_PROVENANCE.md` file.

## Details

The depth-3 asset extends
[`valencia13k_dcst_depth2_merged`](https://pgajer.github.io/linf/reference/valencia13k_dcst_depth2_merged.md)
by one additional hierarchical dCST refinement. All levels use `n0 = 50`
and `low.freq.policy = "absorb"`. This object is intended as a reusable
source for selecting richer VALENCIA-derived component sets without
recomputing the full hierarchy from the 13k source matrix.

## Examples

``` r
data(valencia13k_dcst_depth3_merged)
nrow(valencia13k_dcst_depth3_merged$assignments)
#> [1] 12881
head(valencia13k_dcst_depth3_merged$summaries$depth3)
#>   depth                                                           dcst_label
#> 1     3 Lactobacillus_crispatus__Lactobacillus_iners__Lactobacillus_jensenii
#> 2     3        Gardnerella_vaginalis__Lactobacillus_iners__Atopobium_vaginae
#> 3     3     Lactobacillus_crispatus__Lactobacillus_jensenii__g_Lactobacillus
#> 4     3 Lactobacillus_iners__Lactobacillus_crispatus__Lactobacillus_jensenii
#> 5     3           Atopobium_vaginae__Gardnerella_vaginalis__g_Fastidiosipila
#> 6     3           Gardnerella_vaginalis__Atopobium_vaginae__g_Fastidiosipila
#>     n       prop path_length
#> 1 442 0.03431411           3
#> 2 416 0.03229563           3
#> 3 372 0.02887975           3
#> 4 323 0.02507569           3
#> 5 308 0.02391119           3
#> 6 295 0.02290195           3
```
