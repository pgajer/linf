#' Extended homogeneous-coordinate hypercube embedding
#'
#' @description
#' Computes the zero-aware hypercube embedding associated with one reference
#' component of a nonnegative compositional matrix. For rows with positive
#' reference component, the function forms the ordinary homogeneous ratios
#' against that reference and radially maps them into the unit cube. For rows
#' whose reference component is zero, it uses the L-infinity boundary extension
#' so that the embedding remains defined.
#'
#' @details
#' Let \eqn{x = (x_1,\ldots,x_p)} be a nonnegative row and let \eqn{k} be the
#' reference component. When \eqn{x_k > 0}, define
#' \eqn{z = x_{-k}/x_k}. The embedded row is
#' \deqn{
#'   \sigma_\lambda(\|z\|_1)\frac{z}{\|z\|_\infty},
#'   \qquad
#'   \sigma_\lambda(t) = 1 - \exp(-\lambda t).
#' }
#' When \eqn{x_k = 0}, the embedded row is the L-infinity-normalized boundary
#' vector
#' \deqn{
#'   x_{-k}/\|x_{-k}\|_\infty.
#' }
#' All-zero rows are mapped to all-zero embedded rows by convention.
#'
#' If neither `lambda` nor `log.lambda` is supplied, the scale is chosen from positive-reference
#' rows so that `sigma.target` is attained at the `sigma.quantile` quantile of
#' \eqn{\|z\|_1}. This is a numerical scaling convention for finite datasets; it
#' does not change the reference component or the boundary extension rule.
#' Direction is computed from scaled non-reference values, and radial products
#' and the type-7 quantile are evaluated using logarithms to avoid overflow.
#' The radial map uses `expm1` to preserve very small positive coordinates.
#' Norms beyond the ordinary numeric range are included in automatic calibration.
#'
#' The `log.lambda` attribute always records the natural logarithm of the scale.
#' The `lambda` attribute is `NA_real_` if that scale cannot be represented as a
#' positive finite number. Reuse any fitted scale with
#' `log.lambda = attr(reference.embedding, "log.lambda")`; this also works for
#' extreme calibration data. The all-zero/reference-only fallback scale is one.
#'
#' @param X Nonnegative numeric matrix with samples in rows and features in
#'   columns.
#' @param reference Reference component. May be a column index, feature ID, or
#'   feature label.
#' @param lambda Positive numeric scalar. If `NULL`, choose a data-scaled value
#'   using `sigma.quantile` and `sigma.target`.
#' @param sigma.quantile Quantile of positive finite-reference \eqn{\|z\|_1}
#'   values used for automatic scaling.
#' @param sigma.target Target value of \eqn{\sigma_\lambda(t)} at the selected
#'   quantile during automatic scaling.
#' @param feature.ids Optional stable feature identifiers, length `ncol(X)`.
#' @param feature.labels Optional display labels, length `ncol(X)`.
#' @param tol Nonnegative tolerance. Reference entries `<= tol` are treated as
#'   zero, and L-infinity norms `<= tol` are treated as zero.
#' @param backend Matrix backend: `"auto"`, `"dense"`, or `"sparse"`. Sparse
#'   inputs are accepted, but this implementation returns a dense matrix.
#'   Non-reference zeros remain zero; account for dense storage in large inputs.
#' @param log.lambda Optional finite natural logarithm of the radial scale,
#'   as an alternative to `lambda`. Supply at most one of these two arguments.
#'
#' @return A numeric matrix with `nrow(X)` rows and `ncol(X) - 1` columns. The
#'   columns correspond to the non-reference components. Attributes record the
#'   reference component, lambda choice, and finite/boundary row counts.
#'
#' @examples
#' X <- rbind(
#'   c(A = 2, B = 1, C = 1),
#'   c(A = 0, B = 2, C = 1)
#' )
#' linf.hypercube.embedding(X, reference = "A", lambda = log(2))
#'
#' @export
linf.hypercube.embedding <- function(X,
                                     reference,
                                     lambda = NULL,
                                     sigma.quantile = 0.95,
                                     sigma.target = 0.95,
                                     feature.ids = NULL,
                                     feature.labels = NULL,
                                     tol = 0,
                                     backend = c("auto", "dense", "sparse"),
                                     log.lambda = NULL) {
  prep <- linf.prepare.matrix(X, backend = backend, fun.name = "linf.hypercube.embedding")
  X <- prep$X
  backend <- prep$backend
  linf.validate.matrix(X, backend = backend, fun.name = "linf.hypercube.embedding")

  if (ncol(X) < 2L) stop("linf.hypercube.embedding: X must have at least two features")

  if (backend == "sparse") {
    X <- as.matrix(X)
  }

  if (!is.numeric(tol) || length(tol) != 1L || !is.finite(tol) || tol < 0) {
    stop("linf.hypercube.embedding: tol must be a single non-negative finite number")
  }
  if (!is.numeric(sigma.quantile) || length(sigma.quantile) != 1L ||
      !is.finite(sigma.quantile) || sigma.quantile <= 0 || sigma.quantile > 1) {
    stop("linf.hypercube.embedding: sigma.quantile must be in (0, 1]")
  }
  if (!is.numeric(sigma.target) || length(sigma.target) != 1L ||
      !is.finite(sigma.target) || sigma.target <= 0 || sigma.target >= 1) {
    stop("linf.hypercube.embedding: sigma.target must be in (0, 1)")
  }

  meta <- resolve.linf.feature.meta(
    X,
    feature.ids = feature.ids,
    feature.labels = feature.labels
  )
  feature.ids <- meta$feature.ids
  feature.labels <- meta$feature.labels

  ref.idx <- linf.resolve.reference.index(reference, feature.ids, feature.labels)
  ref.id <- feature.ids[[ref.idx]]
  ref.label <- feature.labels[[ref.idx]]
  other.idx <- setdiff(seq_len(ncol(X)), ref.idx)
  other.ids <- feature.ids[other.idx]
  other.labels <- feature.labels[other.idx]

  denom <- X[, ref.idx]
  others <- X[, other.idx, drop = FALSE]
  finite <- denom > tol

  # Separate direction from radius. No reference ratio or unscaled row sum
  # is needed, so finite inputs cannot overflow either intermediate.
  other.max <- apply(others, 1L, max)
  direction <- matrix(0, nrow(others), ncol(others))
  positive <- other.max > 0
  direction[positive, ] <- others[positive, , drop = FALSE] / other.max[positive]
  log.norm1 <- log.norm.inf <- rep(-Inf, nrow(X))
  log.norm.inf[finite] <- log(other.max[finite]) - log(denom[finite])
  log.norm1[finite] <- log.norm.inf[finite] + log(rowSums(direction[finite, , drop = FALSE]))
  log.tol <- log(tol)

  if (!is.null(lambda) && !is.null(log.lambda)) {
    stop("linf.hypercube.embedding: supply only one of lambda and log.lambda")
  }
  lambda.policy <- "fixed"
  if (!is.null(log.lambda)) {
    if (!is.numeric(log.lambda) || length(log.lambda) != 1L || !is.finite(log.lambda)) {
      stop("linf.hypercube.embedding: log.lambda must be a single finite number")
    }
    lambda.policy <- "log.fixed"
  } else if (!is.null(lambda)) {
    if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
      stop("linf.hypercube.embedding: lambda must be a single positive finite number")
    }
    log.lambda <- log(lambda)
  } else {
    lambda.policy <- "quantile"
    positive.logs <- log.norm1[finite & log.norm1 > log.tol]
    log.lambda <- if (length(positive.logs)) {
      log(-log1p(-sigma.target)) - linf.log.quantile(positive.logs, sigma.quantile)
    } else 0
  }
  # log.lambda is authoritative when the corresponding lambda lies outside the
  # representable positive finite range. Do not expose a misleading zero/Inf.
  lambda <- exp(log.lambda)
  if (!is.finite(lambda) || lambda == 0) lambda <- NA_real_

  out <- matrix(0, nrow = nrow(X), ncol = length(other.idx))
  rownames(out) <- rownames(X)
  colnames(out) <- paste0(other.labels, "_rel_", ref.label)

  if (any(finite)) {
    radius <- -expm1(-exp(log.lambda + log.norm1[finite]))
    finite.direction <- direction[finite, , drop = FALSE]
    finite.direction[log.norm.inf[finite] <= log.tol, ] <- 0
    out[finite, ] <- finite.direction * radius
  }

  boundary <- !finite
  if (any(boundary)) {
    boundary.direction <- direction[boundary, , drop = FALSE]
    boundary.direction[other.max[boundary] <= tol, ] <- 0
    out[boundary, ] <- boundary.direction
  }

  attr(out, "reference.index") <- ref.idx
  attr(out, "reference.id") <- ref.id
  attr(out, "reference.label") <- ref.label
  attr(out, "other.ids") <- other.ids
  attr(out, "other.labels") <- other.labels
  attr(out, "lambda") <- lambda
  attr(out, "log.lambda") <- log.lambda
  attr(out, "lambda.policy") <- lambda.policy
  attr(out, "sigma.quantile") <- sigma.quantile
  attr(out, "sigma.target") <- sigma.target
  attr(out, "finite.reference.count") <- sum(finite)
  attr(out, "zero.reference.count") <- sum(!finite)
  out
}

linf.resolve.reference.index <- function(reference, feature.ids, feature.labels) {
  if (missing(reference) || length(reference) != 1L || is.na(reference)) {
    stop("linf.hypercube.embedding: reference must identify exactly one component")
  }

  if (is.numeric(reference)) {
    if (!is.finite(reference) || reference %% 1 != 0 || reference < 1 || reference > length(feature.ids)) {
      stop("linf.hypercube.embedding: numeric reference is out of range")
    }
    return(as.integer(reference))
  }

  ref <- as.character(reference)
  id.match <- which(feature.ids == ref)
  label.match <- which(feature.labels == ref)
  idx <- unique(c(id.match, label.match))

  if (length(idx) != 1L) {
    stop("linf.hypercube.embedding: reference must match exactly one feature ID or label")
  }
  idx
}


# Type-7 quantile with interpolation in the original scale, calculated in logs.
# Interpolating the logs themselves would change the established calibration.
linf.log.quantile <- function(log.values, probability) {
  values <- sort(log.values)
  h <- 1 + (length(values) - 1) * probability
  lower <- floor(h)
  weight <- h - lower
  if (weight == 0) return(values[[lower]])
  lo <- values[[lower]]
  hi <- values[[lower + 1L]]
  hi + log(weight + (1 - weight) * exp(lo - hi))
}
