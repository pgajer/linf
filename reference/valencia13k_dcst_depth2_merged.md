# Valencia 13k merged depth-2 dCST assignments

A lightweight bundled assignment asset containing merged depth-2
Dominant community state type (dCST) labels for the full Valencia 13k
vaginal microbiome training set.

## Usage

``` r
valencia13k_dcst_depth2_merged
```

## Format

A list with five components:

- assignments:

  Data frame with one row per source sample and columns `sample_id`
  (anonymized bundled-data identifier), `source_row` (row number in the
  filtered Valencia 13k source), `Val_CST`, `Val_subCST`, `dcst_depth1`,
  and `dcst_depth2`.

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

The asset is computed from the Valencia 13k compositional matrix after
L-infinity normalization. dCSTs use `n0 = 50` and the merged
`low.freq.policy = "absorb"` view, so samples from low-support
provisional states are reassigned to retained states rather than stored
as explicit rare buckets.

## Examples

``` r
data(valencia13k_dcst_depth2_merged)
nrow(valencia13k_dcst_depth2_merged$assignments)
#> [1] 12881
head(valencia13k_dcst_depth2_merged$summaries$depth2)
#>   depth                                      dcst_label    n       prop
#> 1     2      Lactobacillus_iners__Gardnerella_vaginalis 1512 0.11738219
#> 2     2 Lactobacillus_crispatus__Lactobacillus_jensenii  972 0.07545998
#> 3     2    Lactobacillus_crispatus__Lactobacillus_iners  868 0.06738607
#> 4     2        Gardnerella_vaginalis__Atopobium_vaginae  819 0.06358202
#> 5     2     Lactobacillus_iners__Lactobacillus_jensenii  746 0.05791476
#> 6     2        Lactobacillus_crispatus__g_Lactobacillus  684 0.05310147
#>   path_length
#> 1           2
#> 2           2
#> 3           2
#> 4           2
#> 5           2
#> 6           2
```
