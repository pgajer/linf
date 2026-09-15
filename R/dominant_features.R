#' Dominant-feature assignment
#'
#' @description
#' Assigns each row to the column achieving its maximum.
#'
#' For each sample (row) of a nonnegative matrix, identifies the dominant
#' feature as the column with the maximum value. Samples with the same dominant
#' feature form a depth-1 dominance sample set. Ties are broken by the first
#' maximum (as in \code{max.col(..., ties.method = "first")}). Rows that are
#' all zero are assigned \code{NA}.
#'
#' Feature IDs default to \code{colnames(S)}; if absent, synthetic IDs
#' \code{"V1", "V2", ..., "Vp"} are generated. Display labels default to the
#' feature IDs unless \code{feature.labels} is supplied. To guarantee a 1-1
#' mapping between columns and both IDs and labels, duplicates are
#' disambiguated via \code{make.unique()}.
#'
#' @param S Numeric matrix (samples x features), typically L-infinity-normalized.
#' @param feature.ids Optional character vector of stable feature identifiers,
#'   length \code{ncol(S)}.
#' @param feature.labels Optional character vector of display labels, length
#'   \code{ncol(S)}.
#' @param tie.method Character. How to resolve ties during dominant-feature
#'   assignment.
#' @param return.value Logical. If `TRUE`, include a `value` vector with row maxima.
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} preserves
#'   sparse input and otherwise uses the dense path.
#'
#' @return A list with components:
#' \itemize{
#'   \item \code{index}: integer index of the dominant column per sample (`NA` for all-zero rows)
#'   \item \code{id}: dominant feature ID per sample (`NA` for all-zero rows)
#'   \item \code{label}: dominant column label per sample (`NA` for all-zero rows)
#'   \item \code{id.levels}: full feature ID set after `make.unique(..., sep = "_")`
#'   \item \code{levels}: full column label set after `make.unique(..., sep = "_")`
#'   \item \code{observed.id.levels}: subset of \code{id.levels} that appear in \code{id}
#'   \item \code{observed.levels}: subset of \code{levels} that appear in \code{label}
#'   \item \code{value}: row maxima (only when \code{return.value = TRUE})
#' }
#'
#' @examples
#' # Basic example with named columns
#' S <- rbind(
#'   a = c(A = 10, B = 5,  C = 0),   # -> A
#'   b = c(A = 0,  B = 0,  C = 0),   # -> NA
#'   c = c(A = 1,  B = 4,  C = 4)    # tie -> first max: B
#' )
#' out <- linf.dominant.features(S)
#' out$index
#' out$label
#' out$levels
#' out$observed.levels
#'
#' # Unnamed columns (synthetic labels V1..Vp), duplicate names disambiguated
#' T <- matrix(c(0,2,  3,1,  0,0), nrow = 3, byrow = TRUE)
#' colnames(T) <- c("X", "X")  # duplicates -> X, X_1
#' linf.dominant.features(T)$levels
#'
#' # With L-infinity normalization in a pipeline
#' M <- normalize.linf(S)
#' linf.dominant.features(M)$label
#'
#' @seealso \code{\link{normalize.linf}}, \code{\link{linf.csts}}
#' @export
linf.dominant.features <- function(S,
                                   feature.ids = NULL,
                                   feature.labels = NULL,
                                   tie.method = c("first", "random", "error"),
                                   return.value = FALSE,
                                   backend = c("auto", "dense", "sparse")) {

  tie.method <- match.arg(tie.method)
  backend <- linf.resolve.backend(S, backend)

  if (backend == "sparse") {
    return(linf.dominant.features.sparse(
      S,
      feature.ids = feature.ids,
      feature.labels = feature.labels,
      tie.method = tie.method,
      return.value = return.value
    ))
  }

  prep <- linf.prepare.matrix(S, backend = "dense", fun.name = "linf.dominant.features")
  X <- prep$X
  linf.validate.matrix(X, backend = "dense", fun.name = "linf.dominant.features")

  meta <- resolve.linf.feature.meta(X, feature.ids = feature.ids, feature.labels = feature.labels)
  id.lev <- meta$feature.ids
  lev <- meta$feature.labels

  idx <- rep(NA_integer_, nrow(X))
  val <- apply(X, 1, max)

  nz <- val > 0
  if (any(nz)) {
    idx[nz] <- apply(X[nz, , drop = FALSE], 1, function(r) {
      m <- max(r)
      j <- which(r == m)
      if (length(j) == 1L) return(j)
      if (tie.method == "first") return(j[1L])
      if (tie.method == "random") return(sample(j, 1L))
      stop("linf.dominant.features: tie encountered and tie.method = 'error'")
    })
  }

  id <- ifelse(is.na(idx), NA_character_, id.lev[idx])
  lbl <- ifelse(is.na(idx), NA_character_, lev[idx])

  observed.id.levels <- id.lev[id.lev %in% id[!is.na(id)]]
  observed.levels <- lev[lev %in% lbl[!is.na(lbl)]]

  out <- list(
    index = idx,
    id = id,
    label = lbl,
    id.levels = id.lev,
    levels = lev,
    observed.id.levels = observed.id.levels,
    observed.levels = observed.levels
  )

  if (return.value) {
    out$value <- val
  }

  out
}
