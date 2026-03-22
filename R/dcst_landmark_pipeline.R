#' Landmark-aware DCST pipeline
#'
#' @description
#' Small wrapper that normalizes a nonnegative matrix with
#' \code{\link{normalize.linf}}, computes depth-1 DCSTs with
#' \code{\link{linf.csts}}, refines once to depth 2 with
#' \code{\link{refine.linf.csts}}, and then computes landmark points for both
#' depths with \code{\link{linf.landmarks}}.
#'
#' The function is intentionally explicit rather than highly abstracted. It
#' keeps the intermediate objects visible and returns them together in one list
#' so downstream workflows can inspect the normalized matrix, the depth-1 DCSTs,
#' the depth-2 DCSTs, and the landmark tables without re-running the pipeline.
#'
#' @param X Numeric matrix (samples x features). Must be finite and nonnegative.
#' @param feature.ids Optional character vector of stable feature identifiers,
#'   length \code{ncol(X)}.
#' @param feature.labels Optional character vector of display labels, length
#'   \code{ncol(X)}.
#' @param n0.depth1 Integer >= 1. Minimum cell size for depth-1 DCSTs.
#' @param n0.depth2 Integer >= 1. Minimum cell size for depth-2 refinement.
#' @param refinement.factor Numeric > 0. Auto-refinement threshold multiplier
#'   passed to \code{\link{refine.linf.csts}}.
#' @param sep Character scalar used to join depth-refined CST path tokens.
#' @param low.freq.policy Character. One of \code{"rare"} or \code{"absorb"}.
#'   Controls the active DCST view while still preserving parallel rare/absorb
#'   views inside the returned CST objects.
#' @param rare.label Character scalar for rare buckets.
#' @param depth1.landmark.types Character vector of landmark types for depth 1.
#' @param depth2.landmark.types Character vector of landmark types for depth 2.
#' @param landmark.view Character. One of \code{"absorb"}, \code{"active"}, or
#'   \code{"rare"}. Defaults to \code{"absorb"} because landmark reporting is
#'   typically requested on the absorb view.
#' @param tie.method Character. Tie handling passed through to
#'   \code{\link{linf.csts}} and \code{\link{linf.landmarks}}.
#' @param verbose Logical. Passed to \code{\link{refine.linf.csts}}.
#'
#' @return A named list with components:
#' \itemize{
#'   \item \code{linf.rel}: L-infinity-normalized matrix.
#'   \item \code{dcst.depth1}: depth-1 \code{"linf.csts"} object.
#'   \item \code{dcst.depth2}: depth-2 \code{"linf.csts"} object.
#'   \item \code{landmarks.depth1}: \code{"linf.landmarks"} object at depth 1.
#'   \item \code{landmarks.depth2}: \code{"linf.landmarks"} object at depth 2.
#'   \item \code{params}: compact record of the pipeline arguments used.
#' }
#'
#' @examples
#' X <- rbind(
#'   s1 = c(10, 8, 1),
#'   s2 = c(9, 7, 2),
#'   s3 = c(8, 2, 7),
#'   s4 = c(7, 1, 8),
#'   s5 = c(1, 10, 2),
#'   s6 = c(2, 9, 1)
#' )
#' ids <- c("asv_1", "asv_2", "asv_3")
#' labels <- c("L. iners 1", "Gard. vaginalis 2", "BVAB1 3")
#'
#' out <- linf.dcst.landmark.pipeline(
#'   X,
#'   feature.ids = ids,
#'   feature.labels = labels,
#'   n0.depth1 = 2,
#'   n0.depth2 = 2,
#'   refinement.factor = 2,
#'   low.freq.policy = "rare",
#'   landmark.view = "absorb",
#'   verbose = FALSE
#' )
#'
#' names(out)
#' out$dcst.depth2$cst.depth
#' out$landmarks.depth2$view
#' @export
linf.dcst.landmark.pipeline <- function(
  X,
  feature.ids = NULL,
  feature.labels = NULL,
  n0.depth1 = 50,
  n0.depth2 = 25,
  refinement.factor = 2,
  sep = "__",
  low.freq.policy = c("rare", "absorb"),
  rare.label = "RARE_DOMINANT",
  depth1.landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep"),
  depth2.landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep"),
  landmark.view = c("absorb", "active", "rare"),
  tie.method = c("first", "random", "error"),
  verbose = FALSE
) {
  low.freq.policy <- match.arg(low.freq.policy)
  landmark.view <- match.arg(landmark.view)
  tie.method <- match.arg(tie.method)

  M <- as.matrix(X)
  storage.mode(M) <- "numeric"

  if (nrow(M) == 0L || ncol(M) == 0L) {
    stop("linf.dcst.landmark.pipeline: matrix has zero rows or columns")
  }
  if (any(!is.finite(M))) {
    stop("linf.dcst.landmark.pipeline: non-finite entries found")
  }
  if (any(M < 0, na.rm = TRUE)) {
    stop("linf.dcst.landmark.pipeline: negative entries found")
  }
  if (!is.numeric(n0.depth1) || length(n0.depth1) != 1L || n0.depth1 < 1 || n0.depth1 %% 1 != 0) {
    stop("linf.dcst.landmark.pipeline: n0.depth1 must be integer >= 1")
  }
  if (!is.numeric(n0.depth2) || length(n0.depth2) != 1L || n0.depth2 < 1 || n0.depth2 %% 1 != 0) {
    stop("linf.dcst.landmark.pipeline: n0.depth2 must be integer >= 1")
  }
  if (!is.numeric(refinement.factor) || length(refinement.factor) != 1L ||
      !is.finite(refinement.factor) || refinement.factor <= 0) {
    stop("linf.dcst.landmark.pipeline: refinement.factor must be a finite numeric > 0")
  }
  if (!is.character(sep) || length(sep) != 1L || !nzchar(sep)) {
    stop("linf.dcst.landmark.pipeline: sep must be a non-empty character scalar")
  }
  if (!is.character(rare.label) || length(rare.label) != 1L || !nzchar(rare.label)) {
    stop("linf.dcst.landmark.pipeline: rare.label must be a non-empty character scalar")
  }

  linf.rel <- normalize.linf(M)

  dcst.depth1 <- linf.csts(
    linf.rel,
    feature.ids = feature.ids,
    feature.labels = feature.labels,
    n0 = n0.depth1,
    low.freq.policy = low.freq.policy,
    rare.label = rare.label,
    tie.method = tie.method
  )

  dcst.depth2 <- refine.linf.csts(
    linf.rel,
    dcst.depth1,
    n0 = n0.depth2,
    refinement.factor = refinement.factor,
    sep = sep,
    low.freq.policy = low.freq.policy,
    rare.label = rare.label,
    verbose = verbose
  )

  landmarks.depth1 <- linf.landmarks(
    linf.rel,
    dcst.depth1,
    depth = 1L,
    view = landmark.view,
    landmark.types = depth1.landmark.types,
    tie.method = tie.method
  )

  landmarks.depth2 <- linf.landmarks(
    linf.rel,
    dcst.depth2,
    depth = 2L,
    view = landmark.view,
    landmark.types = depth2.landmark.types,
    tie.method = tie.method
  )

  list(
    linf.rel = linf.rel,
    dcst.depth1 = dcst.depth1,
    dcst.depth2 = dcst.depth2,
    landmarks.depth1 = landmarks.depth1,
    landmarks.depth2 = landmarks.depth2,
    params = list(
      n0.depth1 = as.integer(n0.depth1),
      n0.depth2 = as.integer(n0.depth2),
      refinement.factor = refinement.factor,
      sep = sep,
      low.freq.policy = low.freq.policy,
      rare.label = rare.label,
      landmark.view = landmark.view,
      depth1.landmark.types = depth1.landmark.types,
      depth2.landmark.types = depth2.landmark.types,
      tie.method = tie.method
    )
  )
}
