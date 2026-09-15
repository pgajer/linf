`%||%` <- function(a, b) if (is.null(a)) b else a

#' Filter ASV count matrix by library size (samples) and prevalence (features)
#'
#' @description
#' Filters an ASV **count** matrix by:
#' (1) removing samples with library size below \code{min.lib};
#' (2) keeping features (taxa) that are "present" in at least
#' \code{ceiling(prev.prop * nrow(S.counts))} samples, where presence is
#' defined by either \code{min.count} (counts) or \code{min.rel} (relative).
#' After feature filtering, zero-total samples are dropped and a row-normalized
#' relative-abundance matrix is returned alongside the filtered counts.
#' Fractional count-like values are accepted, but entries must be finite and
#' nonnegative. Thresholds are validated even when overridden by another rule.
#'
#' @param S.counts Numeric matrix (samples x features) of nonnegative counts.
#' @param min.lib Integer >= 0. Minimum library size (row sum of counts) to keep a sample. Default: 1000.
#' @param prev.prop Numeric in (0,1]. Minimum fraction of samples where a feature must be present.
#'   Default: 0.05.
#' @param min.count Integer (>=1) or NULL. Reads to call "present" (ignored if \code{min.rel} set). Default: 2.
#' @param min.rel Numeric in (0,1) or NULL. Relative abundance to call "present" (overrides \code{min.count}). Default: NULL.
#' @param min.feat.total Integer (>=0) or NULL. Optional minimum total reads across all samples per feature. Default: NULL.
#' @param verbose Logical. Print keep/drop summaries. Default: TRUE.
#'
#' @details
#' Sample filtering is applied on raw counts first. Prevalence is computed on the
#' post-sample-filter matrix using either a count or relative rule. A feature is
#' retained if prevalence \eqn{\ge \lceil \text{prev.prop} \times n_{\text{samples}} \rceil}.
#' After feature filtering, samples with zero remaining counts are dropped and
#' the relative matrix \code{rel} is row-normalized. Zero-total rows have no
#' present features under a relative threshold, and are removed at the final
#' sample-filtering stage. They still count in the prevalence denominator if
#' retained by \code{min.lib = 0}.
#'
#' Empty results retain all documented fields. If no samples survive the initial
#' filter, no features are retained and the returned matrices are 0 by 0.
#' Otherwise, features retained by prevalence can remain as columns of a
#' zero-row result; the returned indices always describe the returned matrices.
#'
#' @return A list with elements:
#' \describe{
#'   \item{counts}{Filtered count matrix (samples x features).}
#'   \item{rel}{Row-normalized relative-abundance matrix.}
#'   \item{kept.sample.idx}{Kept sample indices (original order).}
#'   \item{kept.feature.idx}{Kept feature indices (original order).}
#'   \item{prevalence}{Per-feature prevalence counts for kept features.}
#'   \item{thresholds}{List of thresholds actually used.}
#' }
#'
#' @section Conventions:
#' Dot-delimited names; rows are samples, columns are features. Supply raw counts.
#'
#' @examples
#' set.seed(1)
#' S <- matrix(rpois(100 * 20, lambda = 5), nrow = 100, ncol = 20)
#' res <- filter.asv(S, min.lib = 50, prev.prop = 0.1, min.count = 2)
#' dim(res$counts); dim(res$rel)
#' range(rowSums(res$rel))
#'
#' @export
filter.asv <- function(S.counts,
                       min.lib        = 1000,
                       prev.prop      = 0.05,
                       min.count      = 2,
                       min.rel        = NULL,
                       min.feat.total = NULL,
                       verbose        = TRUE) {

  if (!(is.matrix(S.counts) || is.data.frame(S.counts)) ||
      (is.data.frame(S.counts) && !all(vapply(S.counts, is.numeric, logical(1))))) {
    stop("filter.asv: S.counts must be a numeric matrix or numeric data frame")
  }
  S.counts <- as.matrix(S.counts)
  if (!is.numeric(S.counts) || any(!is.finite(S.counts)) || any(S.counts < 0)) {
    stop("filter.asv: counts must be numeric, finite and nonnegative")
  }
  scalar <- function(x, name, lower, upper = Inf, integer = FALSE,
                     lower.open = FALSE, upper.open = FALSE) {
    if (!is.numeric(x) || length(x) != 1L || !is.finite(x) ||
        (if (lower.open) x <= lower else x < lower) ||
        (if (upper.open) x >= upper else x > upper) ||
        (integer && x %% 1 != 0)) {
      stop("filter.asv: invalid ", name, "; see help('filter.asv') for its scalar range")
    }
  }
  scalar(min.lib, "min.lib", 0, integer = TRUE)
  scalar(prev.prop, "prev.prop", 0, 1, lower.open = TRUE)
  if (!is.null(min.count)) scalar(min.count, "min.count", 1, integer = TRUE)
  if (!is.null(min.rel)) scalar(min.rel, "min.rel", 0, 1, lower.open = TRUE, upper.open = TRUE)
  if (!is.null(min.feat.total)) scalar(min.feat.total, "min.feat.total", 0, integer = TRUE)
  if (!is.logical(verbose) || length(verbose) != 1L || is.na(verbose)) {
    stop("filter.asv: verbose must be TRUE or FALSE")
  }

  # --- sample filter on raw counts ---
  libsize <- rowSums(S.counts)
  if (any(!is.finite(libsize))) stop("filter.asv: library totals exceed the finite numeric range")
  keep.samp <- libsize >= min.lib
  if (verbose) message("Samples kept: ", sum(keep.samp), " / ", nrow(S.counts),
                       "  (min.lib = ", min.lib, ")")
  S.counts <- S.counts[keep.samp, , drop = FALSE]

  if (nrow(S.counts) == 0L || ncol(S.counts) == 0L) {
    warning("No data left after sample filtering.")
  }

  # --- prevalence rule ---
  n <- nrow(S.counts)
  if (!is.null(min.rel)) {
    rs <- rowSums(S.counts)
    present.mat <- matrix(FALSE, nrow(S.counts), ncol(S.counts), dimnames = dimnames(S.counts))
    positive <- rs > 0
    present.mat[positive, ] <- S.counts[positive, , drop = FALSE] / rs[positive] >= min.rel
  } else {
    present.mat <- S.counts >= (min.count %||% 1)
  }
  prev <- colSums(present.mat)
  prev.thld <- ceiling(prev.prop * n)
  keep.feat <- n > 0 & prev >= prev.thld

  if (!is.null(min.feat.total)) {
    keep.feat <- keep.feat & (colSums(S.counts) >= min.feat.total)
  }

  if (verbose) {
    message("Features kept: ", sum(keep.feat), " / ", ncol(S.counts),
            "  (prev.prop = ", prev.prop, ", prev.thld = ", prev.thld, ")")
  }

  S.counts <- S.counts[, keep.feat, drop = FALSE]

  # --- drop zero rows; renormalize ---
  lib2 <- rowSums(S.counts)
  keep.samp2 <- lib2 > 0
  if (verbose && any(!keep.samp2)) {
    message("Dropping ", sum(!keep.samp2), " samples with zero counts after feature filtering.")
  }
  S.counts <- S.counts[keep.samp2, , drop = FALSE]

  rs2 <- rowSums(S.counts)
  S.rel <- S.counts
  if (length(rs2)) S.rel <- S.counts / rs2

  if (nrow(S.rel) > 0) {
    s <- rowSums(S.rel)
    stopifnot(all(is.finite(S.rel)), all(abs(s - 1) < 1e-8))
  }

  list(
    counts = S.counts,
    rel = S.rel,
    kept.sample.idx  = which(keep.samp)[keep.samp2],
    kept.feature.idx = which(keep.feat),
    prevalence = prev[keep.feat],
    thresholds = list(min.lib = min.lib,
                      prev.prop = prev.prop,
                      prev.thld = prev.thld,
                      min.count = min.count,
                      min.rel = min.rel,
                      min.feat.total = min.feat.total)
  )
}
