# linf — L∞ cells & truncated CSTs for compositional data

**linf** provides minimal, dependency-light tools for working with L-infinity representations of compositional data:

- `normalize.linf()` — row-wise L∞ normalization (divide each row by its max; zero rows stay zero)
- `linf.cells()` — dominant-feature assignment per sample (**indices & labels**)
- `linf.csts()` — truncated L∞ CSTs: keep large cells (≥ *n0*), reassign others by restricted argmax
- `filter.asv()` — simple ASV filtering by library size & prevalence (optional helper)
- `asv.to.linf.csts()` — counts → filter → L∞ relatives → truncated CSTs (convenience pipeline)

This package accompanies the Linf paper.

## Installation

```r
# install.packages("devtools")
devtools::install_github("pgajer/linf", build_vignettes = TRUE)
```

or

``` r
install.packages("linf")
```

## Quick start

``` r
library(linf)

set.seed(1)

# toy counts (samples x features)
S.counts <- matrix(rpois(10 * 3, 5), nrow = 10, ncol = 3,
                   dimnames = list(paste0("s",1:10), c("A","B","C")))

# L∞ relatives (nonzero rows have max 1; zeros remain zero)
Z <- normalize.linf(S.counts)
apply(Z, 1, max)
#> returns 1 for nonzero rows, 0 for all-zero rows

# L∞ cells: indices + labels
cells <- linf.cells(Z)
table(cells$label, useNA = "ifany")
cells$levels           # full column-ordered label set (1–1 with columns)
cells$observed.levels  # subset that actually appeared

# Truncated CSTs: keep cells with at least n0 samples
res <- linf.csts(Z, n0 = 4)
res$kept.cells.idx
res$kept.cells.lbl
table(res$cell.label, useNA = "ifany")
```

## From counts to truncated CSTs (pipeline)

``` r
# If you use a prevalence/library-size rule first:
filt <- filter.asv(S.counts, min.lib = 1000, prev.prop = 0.05, min.count = 2)
M <- normalize.linf(filt$counts)
cst <- linf.csts(M, n0 = 50)
```

## Vignette

A short tutorial ships with the package:

``` r
browseVignettes("linf")
# or:
vignette("linf-intro", package = "linf")
```

## Notes & conventions

- Dot-delimited function names (e.g., `filter.asv`, `normalize.linf`)
- `linf.cells()` is invariant to positive row scaling (counts vs
  relatives)
- Ties resolve to the **first** maximum (as in
  `max.col(..., ties.method = "first")`)
- All-zero rows get `NA` for both `index` and `label`

## License

MIT © 2025 Pawel Gajer. See `LICENSE` / `LICENSE.md`.

## Citation

If you use this package, please cite the Linf paper (citation entry will
be added here upon publication).
