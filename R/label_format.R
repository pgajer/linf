#' Format feature labels for L-infinity cells and CSTs
#'
#' @description
#' Builds unique display labels from stable feature IDs and taxonomy strings.
#' This is useful when CST computation should operate on stable internal feature
#' identifiers (for example `asv_4`) while reports and figures should use
#' human-readable labels such as `L. iners 4`.
#'
#' Underscores in taxonomy strings are converted to spaces before abbreviation
#' and aliasing. If multiple features share the same display taxon, an index can
#' be appended either from the global feature ID (for example `asv_4 -> 4`) or
#' by within-taxon order.
#'
#' @param feature.ids Character vector of stable feature identifiers.
#' @param taxonomy Character vector of taxonomy strings, same length as
#'   \code{feature.ids}.
#' @param abbreviations Optional named character vector mapping genus tokens to
#'   display abbreviations, for example \code{c(Lactobacillus = "L.")}.
#' @param aliases Optional named character vector mapping full taxonomy strings
#'   to alternate labels, for example \code{c("Ca. Lachnocurva vaginae" = "BVAB1")}.
#' @param duplicate.index One of \code{"global"}, \code{"within_taxon"}, or
#'   \code{"none"}. Controls how duplicate display taxa are disambiguated.
#' @param fallback.to.id Logical. If \code{TRUE}, missing taxonomy values fall
#'   back to the feature ID.
#'
#' @return Character vector of unique display labels.
#' @examples
#' ids <- c("asv_1", "asv_4", "asv_5", "asv_6")
#' tax <- c(
#'   "Lactobacillus iners",
#'   "Lactobacillus iners",
#'   "Megasphaera lornae",
#'   "Ca_Lachnocurva_vaginae"
#' )
#' abbr <- c(
#'   Lactobacillus = "L.",
#'   Gardnerella = "Gard.",
#'   Megasphaera = "Mega."
#' )
#' aliases <- c(
#'   "Ca. Lachnocurva vaginae" = "BVAB1",
#'   "Ca_Lachnocurva_vaginae" = "BVAB1"
#' )
#' format.linf.feature.labels(ids, tax, abbreviations = abbr, aliases = aliases)
#' @export
format.linf.feature.labels <- function(feature.ids,
                                       taxonomy,
                                       abbreviations = NULL,
                                       aliases = NULL,
                                       duplicate.index = c("global", "within_taxon", "none"),
                                       fallback.to.id = TRUE) {
  duplicate.index <- match.arg(duplicate.index)

  if (length(feature.ids) != length(taxonomy)) {
    stop("format.linf.feature.labels: feature.ids and taxonomy must have the same length")
  }

  ids <- as.character(feature.ids)
  tax <- as.character(taxonomy)

  base <- trimws(gsub("_+", " ", tax))
  missing_tax <- is.na(base) | !nzchar(base)
  if (isTRUE(fallback.to.id)) {
    base[missing_tax] <- ids[missing_tax]
  }

  if (!is.null(aliases)) {
    if (is.null(names(aliases)) || any(!nzchar(names(aliases)))) {
      stop("format.linf.feature.labels: aliases must be a named character vector")
    }
    alias_lookup <- aliases
    original_tax <- as.character(taxonomy)
    for (ii in seq_along(base)) {
      key_orig <- original_tax[[ii]]
      key_norm <- base[[ii]]
      if (!is.na(key_orig) && key_orig %in% names(alias_lookup)) {
        base[[ii]] <- alias_lookup[[key_orig]]
      } else if (!is.na(key_norm) && key_norm %in% names(alias_lookup)) {
        base[[ii]] <- alias_lookup[[key_norm]]
      }
    }
  }

  if (!is.null(abbreviations)) {
    if (is.null(names(abbreviations)) || any(!nzchar(names(abbreviations)))) {
      stop("format.linf.feature.labels: abbreviations must be a named character vector")
    }
    for (ii in seq_along(base)) {
      parts <- strsplit(base[[ii]], "\\s+", perl = TRUE)[[1L]]
      if (!length(parts)) next
      genus <- parts[[1L]]
      if (genus %in% names(abbreviations)) {
        parts[[1L]] <- abbreviations[[genus]]
        base[[ii]] <- paste(parts, collapse = " ")
      }
    }
  }

  counts <- table(base)
  display <- base
  duplicated_taxa <- names(counts[counts > 1L])

  if (length(duplicated_taxa)) {
    if (duplicate.index == "global") {
      extracted <- sub("^.*?(\\d+)$", "\\1", ids)
      has_numeric_suffix <- grepl("\\d+$", ids)
      extracted[!has_numeric_suffix] <- as.character(seq_along(ids))[!has_numeric_suffix]
      for (taxon in duplicated_taxa) {
        idx <- which(base == taxon)
        display[idx] <- paste(base[idx], extracted[idx])
      }
    } else if (duplicate.index == "within_taxon") {
      for (taxon in duplicated_taxa) {
        idx <- which(base == taxon)
        display[idx] <- paste(base[idx], seq_along(idx))
      }
    }
  }

  make.unique(display, sep = "_")
}

resolve.linf.feature.meta <- function(X,
                                      feature.ids = NULL,
                                      feature.labels = NULL) {
  ids <- feature.ids
  if (is.null(ids)) {
    ids <- colnames(X)
  }
  if (is.null(ids)) {
    ids <- paste0("V", seq_len(ncol(X)))
  }
  ids <- make.unique(as.character(ids), sep = "_")

  labels <- feature.labels
  if (is.null(labels)) {
    labels <- ids
  }
  if (length(labels) != length(ids)) {
    stop("resolve.linf.feature.meta: feature.labels must have the same length as feature.ids")
  }
  labels <- make.unique(as.character(labels), sep = "_")

  list(feature.ids = ids, feature.labels = labels)
}
