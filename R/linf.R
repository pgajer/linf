# L-infinity normalization utilities

#' L-infinity normalization (row-wise divide by the per-sample maximum)
#'
#' For each sample (row), divides all entries by that row's maximum so the
#' maximum relative abundance equals 1. Rows that are all zero remain all zero.
#'
#' @param S.counts Numeric matrix (samples x features), nonnegative.
#' @return Numeric matrix of the same size; each nonzero row has max == 1.
#'
#' @examples
#' S <- rbind(
#'   c(10, 5, 0),
#'   c(0, 0, 0),
#'   c(1, 2, 4)
#' )
#' Z <- normalize.linf(S)
#' apply(Z, 1, max)          # 1, 0, 1
#' Z[2, ]                     # stays all zeros
#'
#' @export
normalize.linf <- function(S.counts) {
  X <- as.matrix(S.counts)
  storage.mode(X) <- "numeric"
  if (any(!is.finite(X))) stop("normalize.linf: non-finite entries found")
  if (any(X < 0, na.rm = TRUE)) stop("normalize.linf: negative entries found")
  if (nrow(X) == 0L || ncol(X) == 0L) stop("normalize.linf: matrix has zero rows or columns")
  m <- apply(X, 1, max)
  Z <- sweep(X, 1, m, "/")
  if (any(m == 0)) Z[m == 0, ] <- 0
  Z
}

#' Assign L-infinity cell (dominant feature) per sample — indices & labels
#'
#' For each sample (row) of a nonnegative matrix, identifies the \eqn{L^\infty}
#' cell as the column with the maximum value. Ties are broken by the first
#' maximum (as in \code{max.col(..., ties.method = "first")}). Rows that are all
#' zero are assigned \code{NA}.
#'
#' Column labels are taken from \code{colnames(S)}; if absent, synthetic labels
#' \code{"V1", "V2", ..., "Vp"} are generated. To guarantee a 1–1 mapping between
#' columns and labels, duplicate column names are disambiguated via
#' \code{make.unique()}.
#'
#' @param S Numeric matrix (samples x features), nonnegative. Can be raw counts
#'   or relative abundances (e.g., L1- or L-infinity-normalized). Must be finite
#'   (no \code{NA}/\code{NaN}/\code{Inf}) and have positive dimensions.
#'
#' @return A list with four elements:
#' \describe{
#'   \item{\code{index}}{Integer vector (length \code{nrow(S)}). The 1-based
#'     column index of the dominant feature per row; \code{NA_integer_} for
#'     all-zero rows. Row names (if any) are preserved as names.}
#'   \item{\code{label}}{Character vector (length \code{nrow(S)}). The feature
#'     label corresponding to \code{index}; \code{NA_character_} for all-zero
#'     rows. Row names (if any) are preserved as names.}
#'   \item{\code{levels}}{Character vector (length \code{ncol(S)}). The full,
#'     column-ordered set of feature labels (after \code{make.unique}), in 1–1
#'     correspondence with column indices. Use this for reliable label↔index
#'     mapping.}
#'   \item{\code{observed.levels}}{Character vector of the subset of
#'     \code{levels} that actually appear in \code{label}, in column order
#'     (deterministic). Useful for compact summaries.}
#' }
#'
#' @details
#' - The dominant-feature assignment is invariant under positive row scaling, so
#'   counts or relatives are both acceptable inputs. For consistency with
#'   \emph{L-infinity}-based workflows, consider using \code{\link{normalize.linf}}
#'   before calling this function.
#' - All-zero rows are detected via \code{rowSums(S != 0) == 0} and assigned
#'   \code{NA} in both \code{index} and \code{label}.
#'
#' @examples
#' # Basic example with named columns
#' S <- rbind(
#'   a = c(A = 10, B = 5,  C = 0),   # -> A
#'   b = c(A = 0,  B = 0,  C = 0),   # -> NA
#'   c = c(A = 1,  B = 4,  C = 4)    # tie -> first max: B
#' )
#' out <- linf.cells(S)
#' out$index
#' out$label
#' out$levels
#' out$observed.levels
#'
#' # Unnamed columns (synthetic labels V1..Vp), duplicate names disambiguated
#' T <- matrix(c(0,2,  3,1,  0,0), nrow = 3, byrow = TRUE)
#' colnames(T) <- c("X", "X")  # duplicates -> X, X_1
#' linf.cells(T)$levels
#'
#' # With L-infinity normalization in a pipeline
#' M <- normalize.linf(S)
#' linf.cells(M)$label
#'
#' @seealso \code{\link{normalize.linf}}, \code{\link{linf.csts}}
#' @export
linf.cells <- function(S) {
  X <- as.matrix(S)
  storage.mode(X) <- "numeric"
  if (any(!is.finite(X))) stop("linf.cells: non-finite entries found")
  if (any(X < 0, na.rm = TRUE)) stop("linf.cells: negative entries found")
  if (nrow(X) == 0L || ncol(X) == 0L) stop("linf.cells: matrix has zero rows or columns")

  lev <- colnames(X)
  if (is.null(lev)) lev <- paste0("V", seq_len(ncol(X)))
  lev <- make.unique(lev, sep = "_")  # ensure 1–1 with columns

  zero_row <- rowSums(X != 0) == 0
  idx <- max.col(X, ties.method = "first")
  idx[zero_row] <- NA_integer_

  lbl <- rep(NA_character_, length(idx))
  ok <- !is.na(idx)
  lbl[ok] <- lev[idx[ok]]

  names(idx) <- rownames(X)
  names(lbl) <- rownames(X)

  # observed levels in column order (deterministic)
  obs_lev <- lev[lev %in% lbl]

  list(index = idx, label = lbl, levels = lev, observed.levels = obs_lev)
}

#' Truncated L-infinity CSTs (label-aware)
#'
#' Keep L∞ cells with at least n0 samples (by label); reassign samples from small
#' cells to argmax restricted to the kept set. Returns both indices and labels.
#'
#' @param S  Numeric matrix (samples x features), typically L∞-relatives.
#' @param n0 Integer >= 1. Minimum size for a cell to be kept.
#' @return list with:
#'   - cell.index     : integer vector (final truncated indices)
#'   - cell.label     : character vector (final truncated labels)
#'   - kept.cells.idx : integer vector of kept feature indices
#'   - kept.cells.lbl : character vector of kept feature labels
#'   - raw.index      : integer vector (pre-truncation)
#'   - raw.label      : character vector (pre-truncation)
#'   - size.table     : table of original cell sizes (by label, desc)
#' @export
linf.csts <- function(S, n0 = 50) {
  if (!is.numeric(n0) || length(n0) != 1L || n0 < 1 || n0 %% 1 != 0)
    stop("linf.csts: n0 must be integer >= 1")

  X <- as.matrix(S)
  storage.mode(X) <- "numeric"
  if (any(!is.finite(X))) stop("linf.csts: non-finite entries found")
  if (any(X < 0, na.rm = TRUE)) stop("linf.csts: negative entries found")
  if (nrow(X) == 0L || ncol(X) == 0L) stop("linf.csts: matrix has zero rows or columns")

  # Column labels (synthesize if absent) for stable label semantics
  lev <- colnames(X)
  if (is.null(lev)) lev <- paste0("V", seq_len(ncol(X)))

  raw <- linf.cells(X)  # uses same lev
  # Size table by label (exclude NAs), descending
  tab <- sort(table(raw$label[!is.na(raw$label)]), decreasing = TRUE)

  kept.lbl <- names(tab[tab >= n0])
  kept.idx <- match(kept.lbl, lev)

  # If nothing meets threshold: everything NA
  if (length(kept.idx) == 0L) {
    return(list(
      cell.index     = rep(NA_integer_, nrow(X)),
      cell.label     = rep(NA_character_, nrow(X)),
      kept.cells.idx = integer(0L),
      kept.cells.lbl = character(0L),
      raw.index      = raw$index,
      raw.label      = raw$label,
      size.table     = tab
    ))
  }

  # Start from raw, then reassign small-cell samples
  cell.idx  <- raw$index
  cell.lbl  <- raw$label
  small     <- which(!is.na(raw$label) & !(raw$label %in% kept.lbl))

  if (length(small)) {
    # For each small-row, restricted argmax among kept indices (no column dropping)
    reassigned.idx <- apply(X[small, , drop = FALSE], 1, function(r) {
      kvals <- r[kept.idx]
      kept.idx[which.max(kvals)]
    })
    cell.idx[small] <- reassigned.idx
    cell.lbl[small] <- lev[reassigned.idx]
  }

  # Preserve rownames
  names(cell.idx) <- rownames(X)
  names(cell.lbl) <- rownames(X)

  list(
    cell.index     = cell.idx,
    cell.label     = cell.lbl,
    kept.cells.idx = kept.idx,
    kept.cells.lbl = kept.lbl,
    raw.index      = raw$index,
    raw.label      = raw$label,
    size.table     = tab
  )
}

#' From counts -> filter.asv -> L-infinity relatives -> truncated L-infinity CSTs
#'
#' Convenience pipeline: runs \code{filter.asv()} on counts, computes
#' L-infinity-normalized relatives on the filtered counts, and assigns truncated
#' L-infinity CSTs.
#'
#' @param S.counts Numeric matrix of counts (samples x features).
#' @param ... Arguments passed to \code{filter.asv()} (e.g., \code{min.lib},
#'   \code{prev.prop}, \code{min.count}, \code{min.rel}).
#' @return A list with:
#' \itemize{
#'   \item \code{filter}: the full result from \code{filter.asv()}.
#'   \item \code{linf.rel}: L-infinity-normalized matrix on filtered counts.
#'   \item \code{csts}: result from \code{linf.csts()} on \code{linf.rel}.
#' }
#'
#' @examples
#' set.seed(1)
#' S <- matrix(rpois(100, lambda = 5), nrow = 20, ncol = 5)
#' res <- asv.to.linf.csts(S, min.lib = 10, prev.prop = 0.1, min.count = 1)
#' names(res)
#' apply(res$linf.rel, 1, max)  # should be 1 (or 0 for all-zero rows)
#'
#' @export
asv.to.linf.csts <- function(S.counts, ...) {
  filt <- filter.asv(S.counts, ...)
  M <- normalize.linf(filt$counts)
  cst <- linf.csts(M)
  list(filter = filt, linf.rel = M, csts = cst)
}
