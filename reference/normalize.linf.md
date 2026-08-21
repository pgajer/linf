# L-infinity normalization (row-wise)

Scales each row of a numeric matrix by its L-infinity norm (row
maximum). Rows whose maximum is zero (or below tolerance) are left
unchanged and remain all-zero.

## Usage

``` r
normalize.linf(X, tol = 0, backend = c("auto", "dense", "sparse"))
```

## Arguments

- X:

  Numeric matrix (samples x features).

- tol:

  Numeric \>= 0. Values with row max \<= tol are treated as zero rows.
  Default: 0 (exact zero only).

- backend:

  Character. Matrix backend to use: `"auto"`, `"dense"`, or `"sparse"`.
  The default `"auto"` preserves sparse input and otherwise uses the
  dense path.

## Value

Numeric matrix of same dimensions as X, L-infinity normalized.

## Details

Zero rows have undefined L-infinity direction. By convention, they are
preserved as all-zero rows and will yield NA labels in downstream
dominant-feature or dCST assignment.

## Examples

``` r
X <- rbind(
  sample1 = c(A = 2, B = 1, C = 0),
  sample2 = c(A = 0, B = 0, C = 0)
)
Z <- normalize.linf(X)
Z
#>         A   B C
#> sample1 1 0.5 0
#> sample2 0 0.0 0
apply(Z, 1, max)
#> sample1 sample2 
#>       1       0 
```
