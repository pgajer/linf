# Dominant-feature assignment

Assigns each row to the column achieving its maximum.

For each sample (row) of a nonnegative matrix, identifies the dominant
feature as the column with the maximum value. Samples with the same
dominant feature form a depth-1 dominance sample set. Ties are broken by
the first maximum (as in `max.col(..., ties.method = "first")`). Rows
that are all zero are assigned `NA`.

Feature IDs default to `colnames(S)`; if absent, synthetic IDs
`"V1", "V2", ..., "Vp"` are generated. Display labels default to the
feature IDs unless `feature.labels` is supplied. To guarantee a 1-1
mapping between columns and both IDs and labels, duplicates are
disambiguated via
[`make.unique()`](https://rdrr.io/r/base/make.unique.html).

## Usage

``` r
linf.dominant.features(
  S,
  feature.ids = NULL,
  feature.labels = NULL,
  tie.method = c("first", "random", "error"),
  return.value = FALSE,
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- S:

  Numeric matrix (samples x features), typically L-infinity-normalized.

- feature.ids:

  Optional character vector of stable feature identifiers, length
  `ncol(S)`.

- feature.labels:

  Optional character vector of display labels, length `ncol(S)`.

- tie.method:

  Character. How to resolve ties during dominant-feature assignment.

- return.value:

  Logical. If `TRUE`, include a `value` vector with row maxima.

- backend:

  Character. Matrix backend to use: `"auto"`, `"dense"`, or `"sparse"`.
  The default `"auto"` preserves sparse input and otherwise uses the
  dense path.

## Value

A list with components:

- `index`: integer index of the dominant column per sample (`NA` for
  all-zero rows)

- `id`: dominant feature ID per sample (`NA` for all-zero rows)

- `label`: dominant column label per sample (`NA` for all-zero rows)

- `id.levels`: full feature ID set after `make.unique(..., sep = "_")`

- `levels`: full column label set after `make.unique(..., sep = "_")`

- `observed.id.levels`: subset of `id.levels` that appear in `id`

- `observed.levels`: subset of `levels` that appear in `label`

- `value`: row maxima (only when `return.value = TRUE`)

## See also

[`normalize.linf`](https://pgajer.github.io/linf/reference/normalize.linf.md),
[`linf.csts`](https://pgajer.github.io/linf/reference/linf.csts.md)

## Examples

``` r
# Basic example with named columns
S <- rbind(
  a = c(A = 10, B = 5,  C = 0),   # -> A
  b = c(A = 0,  B = 0,  C = 0),   # -> NA
  c = c(A = 1,  B = 4,  C = 4)    # tie -> first max: B
)
out <- linf.dominant.features(S)
out$index
#> [1]  1 NA  2
out$label
#> [1] "A" NA  "B"
out$levels
#> [1] "A" "B" "C"
out$observed.levels
#> [1] "A" "B"

# Unnamed columns (synthetic labels V1..Vp), duplicate names disambiguated
T <- matrix(c(0,2,  3,1,  0,0), nrow = 3, byrow = TRUE)
colnames(T) <- c("X", "X")  # duplicates -> X, X_1
linf.dominant.features(T)$levels
#> [1] "X"   "X_1"

# With L-infinity normalization in a pipeline
M <- normalize.linf(S)
linf.dominant.features(M)$label
#> [1] "A" NA  "B"
```
