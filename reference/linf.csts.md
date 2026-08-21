# Truncated dominant community state types with configurable low-support handling

Forms provisional depth-1 dominance sample sets from the dominant
feature of each sample and then applies the minimum support threshold
`n0`. Sets with fewer than `n0` samples are handled according to
`low.freq.policy`:

- `"pure"`: retain only sets with support \>= `n0` as named dCSTs and
  collapse all low-support sets into `rare.label`.

- `"absorb"`: reassign each low-support sample to the retained state
  with the largest value among the retained features (ties handled by
  `tie.method`).

## Usage

``` r
linf.csts(
  S,
  feature.ids = NULL,
  feature.labels = NULL,
  n0 = 50,
  low.freq.policy = c("pure", "absorb"),
  rare.label = "RARE_DOMINANT",
  tie.method = c("first", "random", "error"),
  return.diagnostics = FALSE,
  return.landmarks = FALSE,
  landmark.types = c("endpoint.max", "endpoint.min"),
  landmark.view = c("active", "pure", "absorb"),
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- S:

  Numeric matrix (samples x features), typically L-infinity relatives.

- feature.ids:

  Optional character vector of stable feature identifiers, length
  `ncol(S)`.

- feature.labels:

  Optional character vector of display labels, length `ncol(S)`.

- n0:

  Integer \>= 1. Minimum support required to retain a dominance sample
  set.

- low.freq.policy:

  Character. One of `"pure"` or `"absorb"`. Default: `"pure"`.

- rare.label:

  Character scalar used when `low.freq.policy = "pure"`. Default:
  `"RARE_DOMINANT"`.

- tie.method:

  Character. Tie handling passed to
  [`linf.dominant.features()`](https://pgajer.github.io/linf/reference/linf.dominant.features.md)
  and used during absorb reassignment ("first", "random", "error").

- return.diagnostics:

  Logical. If TRUE, return reassignment diagnostics.

- return.landmarks:

  Logical. If TRUE, attach a depth-1 landmark summary computed by
  [`linf.landmarks`](https://pgajer.github.io/linf/reference/linf.landmarks.md).

- landmark.types:

  Character vector of landmark types passed to
  [`linf.landmarks`](https://pgajer.github.io/linf/reference/linf.landmarks.md)
  when `return.landmarks = TRUE`.

- landmark.view:

  Character. Landmark view passed to
  [`linf.landmarks`](https://pgajer.github.io/linf/reference/linf.landmarks.md)
  when `return.landmarks = TRUE`.

- backend:

  Character. Matrix backend to use: `"auto"`, `"dense"`, or `"sparse"`.
  The default `"auto"` preserves sparse input and otherwise uses the
  dense path.

## Value

List with:

- `depth1.feature.index`, `depth1.feature.id`, `depth1.feature.label`:
  active depth-1 assignment

- `lineage.id`, `lineage.label`: active leaf-lineage assignment

- policy-specific variants of the depth-1 and leaf-lineage fields,
  ending in `.pure` or `.absorb`

- `lineage.ids`, `lineage.labels`: active hierarchy, plus
  policy-specific `.pure` and `.absorb` hierarchies

- `depth`: current hierarchy depth

- `retained.feature.indices`, `retained.feature.ids`,
  `retained.feature.labels`

- `provisional.feature.index`, `provisional.feature.id`,
  `provisional.feature.label`

- `feature.ids`, `feature.labels`

- `size.table`, `size.table.id`

- `n0`, `low.freq.policy`, `rare.label`

- `diagnostics` (if `return.diagnostics = TRUE`)

- `landmarks` (if `return.landmarks = TRUE`)

## Examples

``` r
X <- rbind(
  s1 = c(A = 10, B = 2, C = 1),
  s2 = c(A = 9, B = 3, C = 1),
  s3 = c(A = 1, B = 10, C = 2),
  s4 = c(A = 1, B = 9, C = 3),
  s5 = c(A = 1, B = 2, C = 10)
)
fit <- linf.csts(normalize.linf(X), n0 = 2, low.freq.policy = "pure")
table(fit$lineage.label)
#> 
#>             A             B RARE_DOMINANT 
#>             2             2             1 
fit$retained.feature.ids
#> [1] "A" "B"
```
