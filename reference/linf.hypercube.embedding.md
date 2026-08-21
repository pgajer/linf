# Extended homogeneous-coordinate hypercube embedding

Computes the zero-aware hypercube embedding associated with one
reference component of a nonnegative compositional matrix. For rows with
positive reference component, the function forms the ordinary
homogeneous ratios against that reference and radially maps them into
the unit cube. For rows whose reference component is zero, it uses the
L-infinity boundary extension so that the embedding remains defined.

## Usage

``` r
linf.hypercube.embedding(
  X,
  reference,
  lambda = NULL,
  sigma.quantile = 0.95,
  sigma.target = 0.95,
  feature.ids = NULL,
  feature.labels = NULL,
  tol = 0,
  backend = c("auto", "dense", "sparse")
)
```

## Arguments

- X:

  Nonnegative numeric matrix with samples in rows and features in
  columns.

- reference:

  Reference component. May be a column index, feature ID, or feature
  label.

- lambda:

  Positive numeric scalar. If `NULL`, choose a data-scaled value using
  `sigma.quantile` and `sigma.target`.

- sigma.quantile:

  Quantile of positive finite-reference \\\\z\\\_1\\ values used when
  `lambda = NULL`.

- sigma.target:

  Target value of \\\sigma\_\lambda(t)\\ at the selected quantile when
  `lambda = NULL`.

- feature.ids:

  Optional stable feature identifiers, length `ncol(X)`.

- feature.labels:

  Optional display labels, length `ncol(X)`.

- tol:

  Nonnegative tolerance. Reference entries `<= tol` are treated as zero,
  and L-infinity norms `<= tol` are treated as zero.

- backend:

  Matrix backend: `"auto"`, `"dense"`, or `"sparse"`. Sparse inputs are
  accepted, but the returned embedding is a dense matrix because
  homogeneous-coordinate embeddings are generally dense.

## Value

A numeric matrix with `nrow(X)` rows and `ncol(X) - 1` columns. The
columns correspond to the non-reference components. Attributes record
the reference component, lambda choice, and finite/boundary row counts.

## Details

Let \\x = (x_1,\ldots,x_p)\\ be a nonnegative row and let \\k\\ be the
reference component. When \\x_k \> 0\\, define \\z = x\_{-k}/x_k\\. The
embedded row is \$\$ \sigma\_\lambda(\\z\\\_1)\frac{z}{\\z\\\_\infty},
\qquad \sigma\_\lambda(t) = 1 - \exp(-\lambda t). \$\$ When \\x_k = 0\\,
the embedded row is the L-infinity-normalized boundary vector \$\$
x\_{-k}/\\x\_{-k}\\\_\infty. \$\$ All-zero rows are mapped to all-zero
embedded rows by convention.

If `lambda` is not supplied, it is chosen from the positive
finite-reference rows so that `sigma.target` is attained at the
`sigma.quantile` quantile of \\\\z\\\_1\\. This is a numerical scaling
convention for finite datasets; it does not change the reference
component or the boundary extension rule.

## Examples

``` r
X <- rbind(
  c(A = 2, B = 1, C = 1),
  c(A = 0, B = 2, C = 1)
)
linf.hypercube.embedding(X, reference = "A", lambda = log(2))
#>      B_rel_A C_rel_A
#> [1,]     0.5     0.5
#> [2,]     1.0     0.5
#> attr(,"reference.index")
#> [1] 1
#> attr(,"reference.id")
#> [1] "A"
#> attr(,"reference.label")
#> [1] "A"
#> attr(,"other.ids")
#> [1] "B" "C"
#> attr(,"other.labels")
#> [1] "B" "C"
#> attr(,"lambda")
#> [1] 0.6931472
#> attr(,"lambda.policy")
#> [1] "fixed"
#> attr(,"sigma.quantile")
#> [1] 0.95
#> attr(,"sigma.target")
#> [1] 0.95
#> attr(,"finite.reference.count")
#> [1] 1
#> attr(,"zero.reference.count")
#> [1] 1
```
