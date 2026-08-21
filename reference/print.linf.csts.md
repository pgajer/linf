# Print method for `"linf.csts"`

Print method for `"linf.csts"`

## Usage

``` r
# S3 method for class 'linf.csts'
print(x, ...)
```

## Arguments

- x:

  A `"linf.csts"` object.

- ...:

  Unused.

## Value

The input object, invisibly.

## Examples

``` r
M <- rbind(
  s1 = c(A = 1.0, B = 0.2),
  s2 = c(A = 0.8, B = 1.0)
)
fit <- linf.csts(M, n0 = 1)
print(fit)
#> 
#> ================================================================================
#> Dominant Community State Type Hierarchy
#> ================================================================================
#> Total samples:  2 
#> Max depth:      1 
#> Low-freq policy: pure (rare.label: RARE_DOMINANT)
#> --------------------------------------------------------------------------------
#> 
#> Depth 1 
#>   A: 1
#>   B: 1
#> 
```
