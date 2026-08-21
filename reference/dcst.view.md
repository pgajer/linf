# Select a stored dCST policy view

Returns a copy of a `"linf.csts"` object using either the stored
`"pure"` or `"absorb"` hierarchy. This does not recompute dCSTs.

## Usage

``` r
dcst.view(csts, view = c("absorb", "pure"))
```

## Arguments

- csts:

  A `"linf.csts"` object produced by
  [`linf.csts`](https://pgajer.github.io/linf/reference/linf.csts.md)
  and optionally refined by repeated calls to
  [`refine.linf.csts`](https://pgajer.github.io/linf/reference/refine.linf.csts.md).

- view:

  Character. The stored policy view to activate: `"pure"` or `"absorb"`.

## Value

A `"linf.csts"` object using the requested view.

## Examples

``` r
M <- rbind(
  s1 = c(A = 1.0, B = 0.2, C = 0.1),
  s2 = c(A = 0.9, B = 0.3, C = 0.1),
  s3 = c(A = 0.2, B = 1.0, C = 0.1),
  s4 = c(A = 0.2, B = 0.1, C = 1.0)
)
fit <- linf.csts(M, n0 = 2, low.freq.policy = "pure")
table(fit$lineage.label)
#> 
#>             A RARE_DOMINANT 
#>             2             2 
table(dcst.view(fit, view = "absorb")$lineage.label)
#> 
#> A 
#> 4 
```
