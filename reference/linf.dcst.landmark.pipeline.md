# Landmark-aware dCST pipeline

Small wrapper that normalizes a nonnegative matrix with
[`normalize.linf`](https://pgajer.github.io/linf/reference/normalize.linf.md),
computes depth-1 dCSTs with
[`linf.csts`](https://pgajer.github.io/linf/reference/linf.csts.md),
refines once to depth 2 with
[`refine.linf.csts`](https://pgajer.github.io/linf/reference/refine.linf.csts.md),
and then computes landmark points for both depths with
[`linf.landmarks`](https://pgajer.github.io/linf/reference/linf.landmarks.md).

The function is intentionally explicit rather than highly abstracted. It
keeps the intermediate objects visible and returns them together in one
list so downstream workflows can inspect the normalized matrix, the
depth-1 dCSTs, the depth-2 dCSTs, and the landmark tables without
re-running the pipeline.

## Usage

``` r
linf.dcst.landmark.pipeline(
  X,
  feature.ids = NULL,
  feature.labels = NULL,
  n0.depth1 = 50,
  n0.depth2 = 25,
  refinement.factor = 2,
  sep = "__",
  low.freq.policy = c("pure", "absorb"),
  rare.label = "RARE_DOMINANT",
  depth1.landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep"),
  depth2.landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep"),
  landmark.view = c("absorb", "active", "pure"),
  tie.method = c("first", "random", "error"),
  verbose = FALSE,
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- X:

  Numeric matrix (samples x features). Must be finite and nonnegative.

- feature.ids:

  Optional character vector of stable feature identifiers, length
  `ncol(X)`.

- feature.labels:

  Optional character vector of display labels, length `ncol(X)`.

- n0.depth1:

  Integer \>= 1. Minimum support for a depth-1 dominance sample set.

- n0.depth2:

  Integer \>= 1. Minimum support for a depth-2 child lineage.

- refinement.factor:

  Numeric \> 0. Auto-refinement threshold multiplier passed to
  [`refine.linf.csts`](https://pgajer.github.io/linf/reference/refine.linf.csts.md).

- sep:

  Character scalar used to join depth-refined dCST path tokens.

- low.freq.policy:

  Character. One of `"pure"` or `"absorb"`. Controls the active dCST
  view while preserving both pure and absorb views inside the returned
  dCST objects.

- rare.label:

  Character scalar for rare buckets.

- depth1.landmark.types:

  Character vector of landmark types for depth 1.

- depth2.landmark.types:

  Character vector of landmark types for depth 2.

- landmark.view:

  Character. One of `"absorb"`, `"active"`, or `"pure"`. Defaults to
  `"absorb"` because landmark reporting is typically requested on the
  absorb view.

- tie.method:

  Character. Tie handling passed through to
  [`linf.csts`](https://pgajer.github.io/linf/reference/linf.csts.md)
  and
  [`linf.landmarks`](https://pgajer.github.io/linf/reference/linf.landmarks.md).

- verbose:

  Logical. Passed to
  [`refine.linf.csts`](https://pgajer.github.io/linf/reference/refine.linf.csts.md).

- backend:

  Character. Matrix backend to use: `"auto"`, `"dense"`, or `"sparse"`.

## Value

A named list with components:

- `linf.rel`: L-infinity-normalized matrix.

- `dcst.depth1`: depth-1 `"linf.csts"` object.

- `dcst.depth2`: depth-2 `"linf.csts"` object.

- `landmarks.depth1`: `"linf.landmarks"` object at depth 1.

- `landmarks.depth2`: `"linf.landmarks"` object at depth 2.

- `params`: compact record of the pipeline arguments used.

## Examples

``` r
X <- rbind(
  s1 = c(10, 8, 1),
  s2 = c(9, 7, 2),
  s3 = c(8, 2, 7),
  s4 = c(7, 1, 8),
  s5 = c(1, 10, 2),
  s6 = c(2, 9, 1)
)
ids <- c("asv_1", "asv_2", "asv_3")
labels <- c("L. iners 1", "Gard. vaginalis 2", "BVAB1 3")

out <- linf.dcst.landmark.pipeline(
  X,
  feature.ids = ids,
  feature.labels = labels,
  n0.depth1 = 2,
  n0.depth2 = 2,
  refinement.factor = 2,
  low.freq.policy = "pure",
  landmark.view = "absorb",
  verbose = FALSE
)

names(out)
#> [1] "linf.rel"         "dcst.depth1"      "dcst.depth2"      "landmarks.depth1"
#> [5] "landmarks.depth2" "params"          
out$dcst.depth2$depth
#> [1] 2
out$landmarks.depth2$view
#> [1] "absorb"
```
