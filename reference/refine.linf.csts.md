# Refine a dCST hierarchy by one level

Selects leaf dominance-lineages and refines them by dropping the
dominant feature(s) encoded in the parent lineage-ID path and
re-applying
[`linf.csts`](https://pgajer.github.io/linf/reference/linf.csts.md) to
the remaining features. The resulting child labels are appended to the
parent label using `sep`.

Low-support child lineages are handled by `low.freq.policy`. When
`low.freq.policy = "pure"`, rare buckets at depth \>= 2 become
parent-prefixed automatically via the hierarchical
`paste(parent, child, sep = sep)`.

## Usage

``` r
refine.linf.csts(
  M,
  csts,
  lineages.to.refine = NULL,
  n0 = 50,
  refinement.factor = 2,
  sep = "__",
  low.freq.policy = c("pure", "absorb"),
  rare.label = "RARE_DOMINANT",
  verbose = TRUE,
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- M:

  Numeric matrix (samples x features) used for refinement. Columns must
  correspond, in order, to the stable feature IDs stored in `csts`.

- csts:

  A `"linf.csts"` object.

- lineages.to.refine:

  Optional character vector of leaf dominance-lineage IDs to refine.
  When `NULL`, lineages are selected automatically using
  `refinement.factor * n0`.

- n0:

  Integer \>= 1. Minimum support required to retain a child lineage
  (passed to `linf.csts`).

- refinement.factor:

  Numeric \> 0. Auto-refine parent lineages with support \>=
  `refinement.factor * n0`.

- sep:

  Character scalar used to concatenate hierarchical labels.

- low.freq.policy:

  Character. One of `"pure"` or `"absorb"`. Default: `"pure"`.

- rare.label:

  Character scalar for rare buckets when `low.freq.policy = "pure"`.

- verbose:

  Logical. If TRUE, emit progress messages.

- backend:

  Character. Matrix backend to use: `"auto"`, `"dense"`, or `"sparse"`.
  The default `"auto"` inherits the backend from `M` or from `csts` when
  available.

## Value

Updated `"linf.csts"` object with `depth` increased by one and updated
`lineage.label`. Policy-specific views are stored in
`lineage.label.pure` and `lineage.label.absorb`.

## Examples

``` r
M <- rbind(
  s1 = c(A = 1.0, B = 0.8, C = 0.1),
  s2 = c(A = 1.0, B = 0.7, C = 0.2),
  s3 = c(A = 1.0, B = 0.6, C = 0.3),
  s4 = c(A = 0.2, B = 1.0, C = 0.8),
  s5 = c(A = 0.1, B = 1.0, C = 0.7),
  s6 = c(A = 0.3, B = 1.0, C = 0.6)
)
depth1 <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")
depth2 <- refine.linf.csts(
  M,
  depth1,
  lineages.to.refine = "A",
  n0 = 2,
  low.freq.policy = "absorb",
  verbose = FALSE
)
depth2$lineage.labels[[2]]
#>     s1     s2     s3     s4     s5     s6 
#> "A__B" "A__B" "A__B"    "B"    "B"    "B" 
```
