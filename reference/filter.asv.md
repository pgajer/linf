# Filter ASV count matrix by library size (samples) and prevalence (features)

Filters an ASV **count** matrix by: (1) removing samples with library
size below `min.lib`; (2) keeping features (taxa) that are "present" in
at least `ceiling(prev.prop * nrow(S.counts))` samples, where presence
is defined by either `min.count` (counts) or `min.rel` (relative). After
feature filtering, zero-total samples are dropped and a row-normalized
relative-abundance matrix is returned alongside the filtered counts.

## Usage

``` r
filter.asv(
  S.counts,
  min.lib = 1000,
  prev.prop = 0.05,
  min.count = 2,
  min.rel = NULL,
  min.feat.total = NULL,
  verbose = TRUE
)
```

## Arguments

- S.counts:

  Numeric matrix (samples x features) of nonnegative counts.

- min.lib:

  Integer. Minimum library size (row sum of counts) to keep a sample.
  Default: 1000.

- prev.prop:

  Numeric in (0,1\]. Minimum fraction of samples where a feature must be
  present. Default: 0.05.

- min.count:

  Integer (\>=1) or NULL. Reads to call "present" (ignored if `min.rel`
  set). Default: 2.

- min.rel:

  Numeric in (0,1) or NULL. Relative abundance to call "present"
  (overrides `min.count`). Default: NULL.

- min.feat.total:

  Integer (\>=0) or NULL. Optional minimum total reads across all
  samples per feature. Default: NULL.

- verbose:

  Logical. Print keep/drop summaries. Default: TRUE.

## Value

A list with elements:

- counts:

  Filtered count matrix (samples x features).

- rel:

  Row-normalized relative-abundance matrix.

- kept.sample.idx:

  Kept sample indices (original order).

- kept.feature.idx:

  Kept feature indices (original order).

- prevalence:

  Per-feature prevalence counts for kept features.

- thresholds:

  List of thresholds actually used.

## Details

Sample filtering is applied on raw counts first. Prevalence is computed
on the post-sample-filter matrix using either a count or relative rule.
A feature is retained if prevalence \\\ge \lceil \text{prev.prop} \times
n\_{\text{samples}} \rceil\\. After feature filtering, samples with zero
remaining counts are dropped and the relative matrix `rel` is
row-normalized.

## Conventions

Dot-delimited names; rows are samples, columns are features. Supply raw
counts.

## Examples

``` r
set.seed(1)
S <- matrix(rpois(100 * 20, lambda = 5), nrow = 100, ncol = 20)
res <- filter.asv(S, min.lib = 50, prev.prop = 0.1, min.count = 2)
#> Samples kept: 100 / 100  (min.lib = 50)
#> Features kept: 20 / 20  (prev.prop = 0.1, prev.thld = 10)
dim(res$counts); dim(res$rel)
#> [1] 100  20
#> [1] 100  20
range(rowSums(res$rel))
#> [1] 1 1
```
