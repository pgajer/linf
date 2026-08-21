# Landmark points for dCST dominance-lineages

Computes representative landmark points for the dominance-lineages of a
`"linf.csts"` object at a chosen depth and view.

Landmark types are defined with respect to the leaf feature of the dCST
path: the last feature ID in the lineage ID path. Lineages whose leaf
token is `rare.label` are reported but skipped for landmark computation
because they do not correspond to a unique target feature.

## Usage

``` r
linf.landmarks(
  M,
  csts,
  depth = NULL,
  view = c("active", "pure", "absorb"),
  landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep"),
  tie.method = c("first", "random", "error"),
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- M:

  Numeric matrix (samples x features) used to build or refine the dCSTs.

- csts:

  A `"linf.csts"` object.

- depth:

  Integer. dCST depth to inspect. Defaults to the leaf depth
  `csts$depth`.

- view:

  Character. One of `"active"`, `"pure"`, or `"absorb"`.

- landmark.types:

  Character vector containing any of `"endpoint.max"`, `"endpoint.min"`,
  `"mean.rep"`, or `"median.rep"`.

- tie.method:

  Character. Tie handling for landmark selection: `"first"`, `"random"`,
  or `"error"`.

- backend:

  Character. Matrix backend to use: `"auto"`, `"dense"`, or `"sparse"`.
  The default `"auto"` preserves sparse input and otherwise uses the
  dense path.

## Value

A list of class `"linf.landmarks"` with components:

- `depth`, `view`, `sep`, `rare.label`

- `feature.ids`, `feature.labels`

- `lineages`: one row per dominance-lineage with computability metadata

- `landmarks`: one row per computed landmark point

## Examples

``` r
M <- rbind(
  s1 = c(A = 1.0, B = 0.2),
  s2 = c(A = 0.9, B = 0.4),
  s3 = c(A = 0.3, B = 1.0),
  s4 = c(A = 0.1, B = 0.9)
)
fit <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")
landmarks <- linf.landmarks(
  M,
  fit,
  landmark.types = c("endpoint.max", "mean.rep")
)
landmarks$landmarks
#>   lineage.id lineage.label landmark.type point.index point.name
#> 1          A             A  endpoint.max           1         s1
#> 2          A             A      mean.rep           2         s2
#> 3          B             B  endpoint.max           3         s3
#> 4          B             B      mean.rep           4         s4
#>   target.feature.id target.feature.label observed.value target.value
#> 1                 A                    A            1.0         1.00
#> 2                 A                    A            0.9         0.95
#> 3                 B                    B            1.0         1.00
#> 4                 B                    B            0.9         0.95
#>   abs.deviation
#> 1          0.00
#> 2          0.05
#> 3          0.00
#> 4          0.05
```
