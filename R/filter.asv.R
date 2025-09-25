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
#'
#' @param S.counts Numeric matrix (samples x features) of nonnegative counts.
#' @param min.lib Integer. Minimum library size (row sum of counts) to keep a sample. Default: 1000.
#' @param prev.prop Numeric in (0,1]. Minimum fraction of samples where a feature must be present.
#'   Default: 0.05.
#' @param min.count Integer (>=1) or NULL. Reads to call "present" (ignored if \code{min.rel} set). Default: 2.
#' @param min.rel Numeric in (0,1) or NULL. Relative abundance to call "present" (overrides \code{min.count}). Default: NULL.
#' @param min.feat.total Integer (>=0) or NULL. Optional minimum total reads across all samples per feature. Default: NULL.
#' @param verbose Logical. Print keep/drop summaries. Default: TRUE.
#'
#' @details
#' Sample filtering is applied on raw counts first. Prevalence is computed on the
#' post–sample-filter matrix using either a count or relative rule. A feature is
#' retained if prevalence \eqn{\ge \lceil \text{prev.prop} \times n_{\text{samples}} \rceil}.
#' After feature filtering, samples with zero remaining counts are dropped and
#' the relative matrix \code{rel} is row-normalized.
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

  stopifnot(is.matrix(S.counts) || is.data.frame(S.counts))
  S.counts <- as.matrix(S.counts)
  storage.mode(S.counts) <- "numeric"

  # --- sample filter on raw counts ---
  libsize <- rowSums(S.counts)
  keep.samp <- libsize >= min.lib
  if (verbose) message("Samples kept: ", sum(keep.samp), " / ", nrow(S.counts),
                       "  (min.lib = ", min.lib, ")")
  S.counts <- S.counts[keep.samp, , drop = FALSE]

  if (nrow(S.counts) == 0L || ncol(S.counts) == 0L) {
    warning("No data left after sample filtering.")
    return(list(counts = S.counts, rel = S.counts,
                kept.sample.idx = integer(0), kept.feature.idx = integer(0),
                prevalence = numeric(0)))
  }

  # --- prevalence rule ---
  n <- nrow(S.counts)
  if (!is.null(min.rel)) {
    rs <- rowSums(S.counts)
    present.mat <- sweep(S.counts, 1, rs, "/") >= min.rel
  } else {
    present.mat <- S.counts >= (min.count %||% 1)
  }
  prev <- colSums(present.mat)
  prev.thld <- ceiling(prev.prop * n)
  keep.feat <- prev >= prev.thld

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
  S.rel <- sweep(S.counts, 1, rs2, "/")

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
