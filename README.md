
<!-- README.md is generated from README.Rmd. Please edit that file -->

<!-- badges: start -->

[![R-CMD-check](https://github.com/pgajer/linf/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pgajer/linf/actions/workflows/R-CMD-check.yaml)
[![DOI](https://zenodo.org/badge/DOI/10.48550/arXiv.2503.21543.svg)](https://doi.org/10.48550/arXiv.2503.21543)
<!-- badges: end -->

# linf: dominant-feature groups for compositional data

**linf** scales each sample by its largest feature, preserves zeros, and
groups samples by supported dominant and subdominant features. Supply a
nonnegative matrix with **samples in rows and features in columns**.
Stop after normalization, fit a hierarchy, or assign new samples to a
reference you have already fitted.

## A first result

After [installation](#installation), this small example shows what the
support threshold changes. The seed makes the toy counts reproducible.

``` r
library(linf)
set.seed(1)
counts <- matrix(rpois(30, 5), nrow = 10,
                 dimnames = list(paste0("s", 1:10), c("A", "B", "C")))
M <- normalize.linf(counts)
dominant <- linf.dominant.features(M)
fit <- linf.csts(M, n0 = 4, low.freq.policy = "absorb")
```

| Feature | Provisional samples | Assigned after absorption |
|:--------|--------------------:|--------------------------:|
| A       |                   4 |                         6 |
| B       |                   4 |                         4 |
| C       |                   2 |                         0 |

A and B each have four provisional samples and meet `n0 = 4`. C has two,
so its samples are reassigned to the retained feature with greatest
abundance: A finishes with six samples and B with four. Absorption does
not mean A was those samples’ original maximum. `M` alone is sufficient
when you only need normalized profiles; it has row maximum one here and
preserves zeros.

## Installation

The CRAN release is **0.3.0**; the source here is **0.3.1 in
development**. The first example works with either. The two new
navigation/dataset guides and the latest identity corrections are
available in the development version.

``` r
# Released package (two installed workflows)
install.packages("linf")

# Development package (four installed vignettes)
# install.packages("devtools")
devtools::install_github("pgajer/linf", build_vignettes = TRUE)
```

After installing the development package, read the rendered guides
locally:

``` r
help("linf", package = "linf")
vignette("function-guide", package = "linf")
vignette("example-datasets", package = "linf")
```

The [published website](https://pgajer.github.io/linf/) is updated
separately and may lag the development source. Use installed help for
the version you run.

## Choose the next step

- **Prepare profiles:** `filter.asv()` filters counts;
  `normalize.linf()` scales rows above tolerance while leaving zero and
  below-tolerance rows unchanged.
- **Fit groups:** `linf.dominant.features()` names each maximum;
  `linf.csts()` adds a support threshold; `refine.linf.csts()`
  subdivides selected groups.
- **Inspect results:** `summary()` accounts for assigned, rare and
  unassigned samples. `linf.landmarks()` selects observed rows by a
  target feature’s value, including rows closest to its mean or median;
  it does not average profiles.
- **Reuse a fit:** `transfer.dcsts()` assigns query samples to a frozen
  hierarchy.

L-infinity normalization is neither unit-sum normalization nor a
log-ratio transform. It leaves zeros intact. Interpret grouping with the
chosen feature set, support thresholds, low-support policy and tie rule;
deeper is not necessarily better. The [function guide
source](vignettes/function-guide.Rmd) explains those choices and gives
an executable workflow with stopping points.

## Gut Microbiome Demonstration

The figure below illustrates depth-1 dCSTs fitted under the absorb
policy to a bundled set of 766 gut microbiome samples from the American
Gut Project (AGP). After filtering, 763 samples and 307 taxa remain.
With *n*₀ = 30, samples from provisional dominance sample sets below the
support threshold are reassigned to the retained state for which they
have the largest normalized abundance; no composite rare category is
shown.

The subset was deliberately stratified to include every sample assigned
to four selected uncommon dCSTs; the remaining slots are a seed-42
simple random sample from the eligible background. Phenotypes do not
influence selection. The object is suitable for demonstrating the
package workflow, but its phenotype frequencies, effect sizes, and
p-values must not be interpreted as population estimates because
inclusion probabilities differ by dCST.

### dCST Size Distribution

<img src="man/figures/readme-dcst-barplot.png" alt="Barplot of absorb-policy depth-1 dCST sizes in 763 filtered AGP gut samples" width="700" />

The 763 filtered samples are assigned among seven retained dCSTs. The
largest are *Bacteroides* (239 samples), *Escherichia-Shigella* (134),
and *Staphylococcus* (102). In total, 159 samples from low-support
provisional dominance sample sets are absorbed into retained states.

See `vignette("linf-intro", package = "linf")` for a reproducible
demonstration using the bundled data.

## Vignettes

The development package includes four installed vignettes. The first two
links open maintained sources on GitHub; the commands open rendered
installed guides. CRAN 0.3.0 contains the two longer workflows only:

- [Finding your way around linf](vignettes/function-guide.Rmd) —
  task-oriented function catalog, matrix conventions, policies,
  landmarks and frozen-reference transfer. Open locally with
  `vignette("function-guide", package = "linf")`.
- [Example datasets and reproducible
  workflows](vignettes/example-datasets.Rmd) — choose among the five
  bundled objects, align metadata, inspect edge cases and try a separate
  reference/query example. Open with
  `vignette("example-datasets", package = "linf")`.
- [Dominant Community State Types: From Normalization to a Gut
  Demonstration](https://pgajer.github.io/linf/articles/linf-intro.html)
  — normalization, dCST construction and descriptive exploration of the
  stratified AGP gut subset.
- [dCSTs for Vaginal Microbiome
  Data](https://pgajer.github.io/linf/articles/linf-vaginal.html) —
  fitting and descriptive agreement with Valencia labels on samples from
  the classifier’s training source, not held-out validation.

The [interactive embedding
article](https://pgajer.github.io/linf/articles/valencia-hypercube-embedding.html)
is website-only. It is not an installed vignette.

``` r
browseVignettes("linf")
```

## Citation

If you use this package, please cite:

> Gajer, P. & Ravel, J. (2025). A New Approach to Compositional Data
> Analysis using L∞-normalization with Applications to Vaginal
> Microbiome. *arXiv preprint arXiv:2503.21543* \[stat.CO\]. doi:
> [10.48550/arXiv.2503.21543](https://doi.org/10.48550/arXiv.2503.21543)

**BibTeX**

``` bibtex
@article{gajer2025linf,
  title   = {A New Approach to Compositional Data Analysis using
             {$L^{\infty}$}-normalization with Applications to
             Vaginal Microbiome},
  author  = {Gajer, Pawel and Ravel, Jacques},
  year    = {2025},
  eprint  = {2503.21543},
  archivePrefix = {arXiv},
  primaryClass  = {stat.CO},
  journal = {arXiv preprint arXiv:2503.21543},
  doi     = {10.48550/arXiv.2503.21543},
  url     = {https://arxiv.org/abs/2503.21543}
}
```

## License

MIT © 2025 Pawel Gajer. See [LICENSE](LICENSE). Bundled-data sources and
upstream terms are recorded in `inst/DATA_PROVENANCE.md`.
