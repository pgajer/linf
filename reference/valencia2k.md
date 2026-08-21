# Valencia 2k vaginal microbiome dataset

A stratified subsample of 2,000 vaginal samples from the Valencia 13k
CST-classifier training set (France et al. 2020), bundled as an example
dataset for demonstrating dCST construction after L-infinity
normalization.

## Usage

``` r
valencia2k
```

## Format

A list with five components:

- rel:

  Numeric matrix (2000 x 178). Compositional relative abundances; each
  row sums to 1. Rows are samples, columns are taxonomic features.

- cst:

  Data frame (2000 x 3) with columns: `sample_id` (character), `Val_CST`
  (Valencia CST assignment: I, II, III, IV-A, IV-B, IV-C, V),
  `Val_subCST` (Valencia sub-CST assignment: I-A, I-B, II, III-A, III-B,
  IV-A, IV-B, IV-C0, IV-C1, IV-C2, IV-C3, IV-C4, V).

- reads:

  Integer vector of length 2000. Per-sample read counts after taxonomic
  filtering. Use `sweep(valencia2k$rel, 1, valencia2k$reads, "*")` to
  reconstruct a count-like matrix.

- taxa:

  Character vector of 178 taxon names (column names of `rel`).

- source:

  Character string documenting provenance.

## Source

Subsampled from the VALENCIA training data at
<https://github.com/ravel-lab/VALENCIA>. See
`data-raw/build_valencia2k.R` and the installed `DATA_PROVENANCE.md`
file for construction and licensing details.

## Details

The subsample preserves proportional representation of all 13 Valencia
sub-CSTs and was drawn with `set.seed(42)`.

The Valencia CST classifier (France et al. 2020) assigns vaginal
microbiome samples to community state types (CSTs) based on
nearest-centroid classification in relative-abundance space. The
original training set contains 12,881 samples and 178 taxonomic features
after filtering.

This 2,000-sample subsample is intended for vignette demonstrations. The
compositional matrix can be used directly with
[`normalize.linf`](https://pgajer.github.io/linf/reference/normalize.linf.md)
and downstream dCST functions. For workflows that require count-like
input, reconstruct approximate counts using the `reads` vector.

## References

France, M. T., Ma, B., Gajer, P., Brown, S., Humphrys, M. S., Holm, J.
B., Waetjen, L. E., Brotman, R. M., & Ravel, J. (2020). VALENCIA: a
nearest centroid classification method for vaginal microbial communities
based on composition. *Microbiome*, 8(1), 166.
[doi:10.1186/s40168-020-00934-6](https://doi.org/10.1186/s40168-020-00934-6)

## Examples

``` r
data(valencia2k)
dim(valencia2k$rel)           # 2000 x 178
#> [1] 2000  178
table(valencia2k$cst$Val_CST) # CST distribution
#> 
#>    I   II  III IV-A IV-B IV-C    V 
#>  557   69  582  141  454  117   80 
head(valencia2k$taxa, 10)     # first 10 taxon names
#>  [1] "Lactobacillus_iners"     "Lactobacillus_crispatus"
#>  [3] "Gardnerella_vaginalis"   "Lactobacillus_jensenii" 
#>  [5] "Lactobacillus_gasseri"   "Atopobium_vaginae"      
#>  [7] "BVAB1"                   "Sneathia_sanguinegens"  
#>  [9] "g_Megasphaera"           "g_Streptococcus"        
```
