# Format feature labels for dominant-feature assignments and dCSTs

Builds unique display labels from stable feature IDs and taxonomy
strings. This is useful when dCST computation should operate on stable
internal feature identifiers (for example `asv_4`) while reports and
figures should use human-readable labels such as `L. iners 4`.

Underscores in taxonomy strings are converted to spaces before
abbreviation and aliasing. If multiple features share the same display
taxon, an index can be appended either from the global feature ID (for
example `asv_4 -> 4`) or by within-taxon order.

## Usage

``` r
linf.feature.labels(
  feature.ids,
  taxonomy,
  abbreviations = NULL,
  aliases = NULL,
  duplicate.index = c("global", "within_taxon", "none"),
  fallback.to.id = TRUE
)
```

## Arguments

- feature.ids:

  Character vector of stable feature identifiers.

- taxonomy:

  Character vector of taxonomy strings, same length as `feature.ids`.

- abbreviations:

  Optional named character vector mapping genus tokens to display
  abbreviations, for example `c(Lactobacillus = "L.")`.

- aliases:

  Optional named character vector mapping full taxonomy strings to
  alternate labels, for example
  `c("Ca. Lachnocurva vaginae" = "BVAB1")`.

- duplicate.index:

  One of `"global"`, `"within_taxon"`, or `"none"`. Controls how
  duplicate display taxa are disambiguated.

- fallback.to.id:

  Logical. If `TRUE`, missing taxonomy values fall back to the feature
  ID.

## Value

Character vector of unique display labels.

## Examples

``` r
ids <- c("asv_1", "asv_4", "asv_5", "asv_6")
tax <- c(
  "Lactobacillus iners",
  "Lactobacillus iners",
  "Megasphaera lornae",
  "Ca_Lachnocurva_vaginae"
)
abbr <- c(
  Lactobacillus = "L.",
  Gardnerella = "Gard.",
  Megasphaera = "Mega."
)
aliases <- c(
  "Ca. Lachnocurva vaginae" = "BVAB1",
  "Ca_Lachnocurva_vaginae" = "BVAB1"
)
linf.feature.labels(ids, tax, abbreviations = abbr, aliases = aliases)
#> [1] "L. iners 1"   "L. iners 4"   "Mega. lornae" "BVAB1"       
```
