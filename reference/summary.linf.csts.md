# Summarize dCST hierarchy statistics

Summarize dCST hierarchy statistics

## Usage

``` r
# S3 method for class 'linf.csts'
summary(object, ...)
```

## Arguments

- object:

  A `"linf.csts"` object.

- ...:

  Unused.

## Value

Data frame with depth-wise statistics

## Examples

``` r
M <- rbind(
  s1 = c(A = 1.0, B = 0.2),
  s2 = c(A = 0.9, B = 0.3),
  s3 = c(A = 0.2, B = 1.0)
)
fit <- linf.csts(M, n0 = 1)
summary(fit)
#>   depth n.lineages total.samples mean.size median.size min.size max.size
#> 1     1          2             3       1.5         1.5        1        2
```
