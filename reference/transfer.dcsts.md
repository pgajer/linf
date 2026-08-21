# Transfer samples into a frozen dCST hierarchy

Assigns rows of a new sample-by-feature matrix into the realized lineage
sets of a fitted `"linf.csts"` hierarchy. The hierarchy is not refit:
each sample is walked through the frozen tree by choosing, at each
depth, the realized child lineage whose newly added feature has the
largest abundance in that sample.

## Usage

``` r
transfer.dcsts(
  X,
  csts,
  depth = NULL,
  view = c("absorb", "active", "pure"),
  match.by = c("feature.ids", "feature.labels"),
  feature.ids = NULL,
  feature.labels = NULL,
  tie.method = c("support", "first", "random", "error"),
  carry.forward.terminal.depths = TRUE,
  sep = NULL,
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- X:

  Nonnegative sample-by-feature matrix to transfer. Columns may be in
  any order and may include features not present in `csts`; missing
  retained features are treated as zero.

- csts:

  A fitted `"linf.csts"` object produced by
  [`linf.csts`](https://pgajer.github.io/linf/reference/linf.csts.md)
  and optionally refined by repeated calls to
  [`refine.linf.csts`](https://pgajer.github.io/linf/reference/refine.linf.csts.md).

- depth:

  Integer vector of requested depths. Defaults to all fitted depths in
  `csts`.

- view:

  Which fitted hierarchy view to transfer into. `"absorb"` uses
  `lineage.labels.absorb`; `"pure"` uses `lineage.labels.pure`;
  `"active"` uses `lineage.labels`.

- match.by:

  Whether columns in `X` are aligned to `csts$feature.ids` or
  `csts$feature.labels`.

- feature.ids:

  Optional feature identifiers for columns of `X`. Defaults to
  `colnames(X)`.

- feature.labels:

  Optional feature labels for columns of `X`. Defaults to `feature.ids`.

- tie.method:

  How to resolve ties among frozen child lineages with equal sample
  abundance. `"support"` chooses the tied lineage with largest reference
  support and then lexical order; `"first"` uses frozen child order;
  `"random"` samples one tied lineage; `"error"` stops.

- carry.forward.terminal.depths:

  Logical. If `TRUE`, an unmatched depth-1 sample may fall back to its
  dominant retained feature. Terminal lineages at later depths are
  assigned only when they appear as realized stable lineage sets in the
  fitted hierarchy.

- sep:

  Separator used in lineage labels. Defaults to `csts$sep` or `"__"`.

- backend:

  Matrix backend; passed to the internal matrix preparer.

## Value

A list of class `"linf.dcst.transfer"` with components:

- `assignment`: character matrix of transferred labels for the requested
  depths.

- `all.depths`: character matrix for all fitted depths.

- `depth`: requested depth vector.

- `view`, `match.by`, `tie.method`: settings used.

## Examples

``` r
X <- rbind(
  s1 = c(A = 10, B = 2, C = 1),
  s2 = c(A = 9, B = 3, C = 1),
  s3 = c(A = 1, B = 10, C = 2),
  s4 = c(A = 1, B = 9, C = 3)
)
M <- normalize.linf(X)
fit <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")
transfer.dcsts(X, fit)$assignment
#>    depth1
#> s1 "A"   
#> s2 "A"   
#> s3 "B"   
#> s4 "B"   
```
