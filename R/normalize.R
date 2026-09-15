# L-infinity normalization utilities

#' L-infinity normalization (row-wise)
#'
#' Scales each row of a numeric matrix by its L-infinity norm (row maximum).
#' Rows whose maximum is at or below tolerance are left unchanged; nonzero
#' entries in these rows are not replaced by zeros.
#'
#' @param X Numeric matrix (samples x features).
#' @param tol Finite numeric >= 0. Rows with maximum <= tol are not scaled.
#'   Default: 0 (exact zero only).
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} preserves
#'   sparse input and otherwise uses the dense path.
#'
#' @return Numeric matrix of same dimensions as X, L-infinity normalized.
#'
#' @details
#' Zero rows have undefined L-infinity direction. By convention, they are
#' preserved as all-zero rows and yield \code{NA} in dominant-feature assignment.
#' The pure dCST view places them in the rare category; the absorb view assigns
#' them to its fallback retained state, if any. A nonzero row left unscaled
#' because of \code{tol} still has a dominant feature.
#'
#' @examples
#' X <- rbind(
#'   sample1 = c(A = 2, B = 1, C = 0),
#'   sample2 = c(A = 0, B = 0, C = 0)
#' )
#' Z <- normalize.linf(X)
#' Z
#' apply(Z, 1, max)
#'
#' @export
normalize.linf <- function(X,
                           tol = 0,
                           backend = c("auto", "dense", "sparse")) {
  prep <- linf.prepare.matrix(X, backend = backend, fun.name = "normalize.linf")
  X <- prep$X
  backend <- prep$backend

  linf.validate.matrix(X, backend = backend, fun.name = "normalize.linf")
  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("normalize.linf: tol must be a single non-negative finite number")
  }

  if (backend == "sparse") {
    return(linf.normalize.sparse(X, tol = tol))
  }

  m <- apply(X, 1, max)

  ## Identify rows with meaningful L-infinity norm
  keep <- m > tol

  if (any(keep)) {
    X[keep, ] <- X[keep, , drop = FALSE] / m[keep]
  }

  X
}
