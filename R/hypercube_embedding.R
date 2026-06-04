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
#' If `lambda` is not supplied, it is chosen from the positive finite-reference
#' rows so that `sigma.target` is attained at the `sigma.quantile` quantile of
#' \eqn{\|z\|_1}. This is a numerical scaling convention for finite datasets; it
#' does not change the reference component or the boundary extension rule.
#'
#' @param X Nonnegative numeric matrix with samples in rows and features in
#'   columns.
#' @param reference Reference component. May be a column index, feature ID, or
#'   feature label.
#' @param lambda Positive numeric scalar. If `NULL`, choose a data-scaled value
#'   using `sigma.quantile` and `sigma.target`.
#' @param sigma.quantile Quantile of positive finite-reference \eqn{\|z\|_1}
#'   values used when `lambda = NULL`.
#' @param sigma.target Target value of \eqn{\sigma_\lambda(t)} at the selected
#'   quantile when `lambda = NULL`.
#' @param feature.ids Optional stable feature identifiers, length `ncol(X)`.
#' @param feature.labels Optional display labels, length `ncol(X)`.
#' @param tol Nonnegative tolerance. Reference entries `<= tol` are treated as
#'   zero, and L-infinity norms `<= tol` are treated as zero.
#' @param backend Matrix backend: `"auto"`, `"dense"`, or `"sparse"`. Sparse
#'   inputs are accepted, but the returned embedding is a dense matrix because
#'   homogeneous-coordinate embeddings are generally dense.
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
                                     backend = c("auto", "dense", "sparse")) {
  prep <- linf.prepare.matrix(X, backend = backend, fun.name = "linf.hypercube.embedding")
  X <- prep$X
  backend <- prep$backend
  linf.validate.matrix(X, backend = backend, fun.name = "linf.hypercube.embedding")

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

  z.norm1 <- rep(NA_real_, nrow(X))
  if (any(finite)) {
    z.norm1[finite] <- rowSums(others[finite, , drop = FALSE] / denom[finite])
  }

  lambda.policy <- "fixed"
  if (is.null(lambda)) {
    lambda.policy <- "quantile"
    positive.norms <- z.norm1[is.finite(z.norm1) & z.norm1 > tol]
    if (length(positive.norms)) {
      q <- stats::quantile(
        positive.norms,
        probs = sigma.quantile,
        names = FALSE,
        na.rm = TRUE,
        type = 7
      )
      lambda <- if (is.finite(q) && q > tol) -log(1 - sigma.target) / q else 1
    } else {
      lambda <- 1
    }
  }
  if (!is.numeric(lambda) || length(lambda) != 1L || !is.finite(lambda) || lambda <= 0) {
    stop("linf.hypercube.embedding: lambda must be a single positive finite number")
  }

  out <- matrix(0, nrow = nrow(X), ncol = length(other.idx))
  rownames(out) <- rownames(X)
  colnames(out) <- paste0(other.labels, "_rel_", ref.label)

  if (any(finite)) {
    z <- others[finite, , drop = FALSE] / denom[finite]
    z.inf <- apply(z, 1L, max)
    z.one <- rowSums(z)
    direction <- matrix(0, nrow = nrow(z), ncol = ncol(z))
    nz <- z.inf > tol
    direction[nz, ] <- z[nz, , drop = FALSE] / z.inf[nz]
    out[finite, ] <- direction * (1 - exp(-lambda * z.one))
  }

  boundary <- !finite
  if (any(boundary)) {
    boundary.x <- others[boundary, , drop = FALSE]
    boundary.inf <- apply(boundary.x, 1L, max)
    boundary.out <- matrix(0, nrow = nrow(boundary.x), ncol = ncol(boundary.x))
    nz <- boundary.inf > tol
    boundary.out[nz, ] <- boundary.x[nz, , drop = FALSE] / boundary.inf[nz]
    out[boundary, ] <- boundary.out
  }

  attr(out, "reference.index") <- ref.idx
  attr(out, "reference.id") <- ref.id
  attr(out, "reference.label") <- ref.label
  attr(out, "other.ids") <- other.ids
  attr(out, "other.labels") <- other.labels
  attr(out, "lambda") <- lambda
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
    if (reference %% 1 != 0 || reference < 1 || reference > length(feature.ids)) {
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
