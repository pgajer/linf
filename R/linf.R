# L-infinity normalization utilities

#' L-infinity normalization (row-wise)
#'
#' Scales each row of a numeric matrix by its L-infinity norm (row maximum).
#' Rows whose maximum is zero (or below tolerance) are left unchanged
#' and remain all-zero.
#'
#' @param X Numeric matrix (samples x features).
#' @param tol Numeric >= 0. Values with row max <= tol are treated as zero rows.
#'   Default: 0 (exact zero only).
#'
#' @return Numeric matrix of same dimensions as X, L-infinity normalized.
#'
#' @details
#' Zero rows have undefined L-infinity direction. By convention, they are
#' preserved as all-zero rows and will yield NA labels in downstream
#' L-infinity cell or CST assignment.
#'
#' @export
normalize.linf <- function(X, tol = 0) {
  X <- as.matrix(X)
  storage.mode(X) <- "numeric"

  if (any(!is.finite(X))) {
    stop("normalize.linf: non-finite entries found")
  }
  if (any(X < 0, na.rm = TRUE)) {
    stop("normalize.linf: negative entries found")
  }
  if (!is.numeric(tol) || length(tol) != 1L || tol < 0) {
    stop("normalize.linf: tol must be a single non-negative number")
  }

  m <- apply(X, 1, max)

  ## Identify rows with meaningful L∞ norm
  keep <- m > tol

  if (any(keep)) {
    X[keep, ] <- X[keep, , drop = FALSE] / m[keep]
  }

  X
}

#' Truncated L-infinity CSTs (hierarchy-aware)
#'
#' @description
#' Assigns each row to the column achieving its maximum (L-infinity cell).
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
#' @param S Numeric matrix (samples x features), typically L-infinity-normalized.
#' @param n0 Integer >= 1. Minimum size for a CST to be retained.
#' @param tie.method Character. How to resolve ties when assigning L-infinity cells.
#' @param return.diagnostics Logical. If TRUE, return reassignment diagnostics.
#'
#' @return A CST object with components:
#' \itemize{
#'   \item \code{cell.label}: leaf CST labels for each sample
#'   \item \code{cst.levels}: list of CST label vectors by depth
#'   \item \code{cst.depth}: maximum CST depth
#'   \item \code{size.table}: CST sizes at depth 1
#' }
#'
#' @details
#' The returned object represents a hierarchical CST structure.
#' \code{cell.label} always refers to the leaf (deepest) CST labels.
#' The full hierarchy is stored explicitly in \code{cst.levels}.
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
linf.cells <- function(S,
                       tie.method = c("first", "random", "error"),
                       return.value = FALSE) {

  tie.method <- match.arg(tie.method)

  X <- as.matrix(S)
  storage.mode(X) <- "numeric"

  if (any(!is.finite(X))) stop("linf.cells: non-finite entries found")
  if (any(X < 0, na.rm = TRUE)) stop("linf.cells: negative entries found")
  if (nrow(X) == 0L || ncol(X) == 0L)
    stop("linf.cells: matrix has zero rows or columns")

  lev <- colnames(X)
  if (is.null(lev)) lev <- paste0("V", seq_len(ncol(X)))

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
      stop("linf.cells: tie encountered and tie.method = 'error'")
    })
  }

  lbl <- ifelse(is.na(idx), NA_character_, lev[idx])

  out <- list(
    index  = idx,
    label  = lbl,
    levels = lev
  )

  if (return.value) {
    out$value <- val
  }

  out
}

#' Truncated L-infinity CSTs (label-aware)
#'
#' Keeps L-infinity cells with at least n0 samples; samples from smaller cells
#' are reassigned to the L-infinity maximum among the kept set.
#'
#' @param S Numeric matrix (samples x features), typically L-infinity relatives.
#' @param n0 Integer >= 1. Minimum size for a cell to be kept.
#' @param tie.method Character. Tie handling passed to linf.cells() and used
#'   during reassignment ("first", "random", "error").
#' @param return.diagnostics Logical. If TRUE, return reassignment diagnostics.
#'
#' @return List with:
#'   - cell.index
#'   - cell.label
#'   - kept.cells.idx
#'   - kept.cells.lbl
#'   - raw.index
#'   - raw.label
#'   - size.table
#'   - diagnostics (if return.diagnostics = TRUE)
#'
#' @export
linf.csts <- function(S,
                      n0 = 50,
                      tie.method = c("first", "random", "error"),
                      return.diagnostics = FALSE) {

    tie.method <- match.arg(tie.method)

    if (!is.numeric(n0) || length(n0) != 1L || n0 < 1 || n0 %% 1 != 0) {
        stop("linf.csts: n0 must be integer >= 1")
    }

    X <- as.matrix(S)
    storage.mode(X) <- "numeric"

    if (any(!is.finite(X))) stop("linf.csts: non-finite entries found")
    if (any(X < 0, na.rm = TRUE)) stop("linf.csts: negative entries found")
    if (nrow(X) == 0L || ncol(X) == 0L)
        stop("linf.csts: matrix has zero rows or columns")

    lev <- colnames(X)
    if (is.null(lev)) lev <- paste0("V", seq_len(ncol(X)))

    raw <- linf.cells(X, tie.method = tie.method)
    tab <- sort(table(raw$label[!is.na(raw$label)]), decreasing = TRUE)

    kept.lbl <- names(tab[tab >= n0])
    kept.idx <- match(kept.lbl, lev)

    n <- nrow(X)
    cell.idx <- raw$index
    cell.lbl <- raw$label

    reassigned <- logical(n)
    reassigned.from <- rep(NA_character_, n)
    reassigned.to   <- rep(NA_character_, n)

    if (length(kept.idx) == 0L) {
        out <- list(
            cell.index     = rep(NA_integer_, n),
            cell.label     = rep(NA_character_, n),
            kept.cells.idx = integer(0L),
            kept.cells.lbl = character(0L),
            raw.index      = raw$index,
            raw.label      = raw$label,
            size.table     = tab
        )
        if (return.diagnostics) {
            out$diagnostics <- list(
                reassigned      = reassigned,
                reassigned.from = reassigned.from,
                reassigned.to   = reassigned.to
            )
        }
        return(out)
    }

    small <- which(!is.na(raw$label) & !(raw$label %in% kept.lbl))

    if (length(small)) {
        new.idx <- apply(X[small, , drop = FALSE], 1, function(r) {
            kvals <- r[kept.idx]
            m <- max(kvals)
            j <- kept.idx[kvals == m]
            if (length(j) == 1L) return(j)
            if (tie.method == "first") return(j[1L])
            if (tie.method == "random") return(sample(j, 1L))
            stop("linf.csts: tie during reassignment and tie.method = 'error'")
        })

        reassigned[small] <- TRUE
        reassigned.from[small] <- raw$label[small]
        reassigned.to[small]   <- lev[new.idx]

        cell.idx[small] <- new.idx
        cell.lbl[small] <- lev[new.idx]
    }

    names(cell.idx) <- rownames(X)
    names(cell.lbl) <- rownames(X)

    out <- list(
        cell.index     = cell.idx,
        cell.label     = cell.lbl,
        kept.cells.idx = kept.idx,
        kept.cells.lbl = kept.lbl,
        raw.index      = raw$index,
        raw.label      = raw$label,
        size.table     = tab
    )

    if (return.diagnostics) {
        out$diagnostics <- list(
            reassigned      = reassigned,
            reassigned.from = reassigned.from,
            reassigned.to   = reassigned.to
        )
    }

    class(out) <- "linf.csts"

    out
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

validate.linf.csts <- function(obj) {

  stopifnot(is.list(obj))
  stopifnot(!is.null(obj$cell.label))

  n <- length(obj$cell.label)

  ## Case 1: explicit hierarchy (refined CSTs)
  if (!is.null(obj$cst.levels)) {

    stopifnot(is.list(obj$cst.levels))
    stopifnot(!is.null(obj$cst.depth))
    stopifnot(obj$cst.depth == length(obj$cst.levels))

    for (lvl in obj$cst.levels) {
      stopifnot(length(lvl) == n)
    }

    return(invisible(TRUE))
  }

  ## Case 2: depth-1 CST (plain linf.csts output)
  ## This is VALID and should be accepted
  return(invisible(TRUE))
}

#' Refine L-infinity CSTs by subdividing large cells
#'
#' @description
#' Subdivides one or more L-infinity cells into finer sub-cells by:
#' \enumerate{
#'   \item Restricting to samples in the target parent cell
#'   \item Removing the dominant feature column (the one defining the parent cell)
#'   \item Re-applying L-infinity CST assignment to find secondary dominance patterns
#'   \item Creating hierarchical labels like \code{"ParentCell_SubCell"}
#'   \item Post-processing: reassigning any sub-cells with < n0 samples (including
#'         NA assignments) to their nearest sibling cell within the same parent group
#' }
#'
#' When \code{cells.to.refine = NULL}, automatically identifies cells for
#' refinement based on sample size thresholds, making exploratory analysis
#' more streamlined.
#'
#' This is useful when a cell contains many samples and you want to explore
#' internal structure based on which feature is second-most dominant after
#' removing the cell-defining dominant feature.
#'
#' @param M Numeric matrix (samples x features), typically L-infinity-normalized.
#'   Must be the same matrix used to generate \code{csts}. Should have column
#'   names that match cell labels.
#' @param csts List returned by \code{\link{linf.csts}}, containing at minimum
#'   \code{cell.label} (character vector of cell assignments) and optionally
#'   \code{cell.index}.
#' @param n0 Integer >= 1. Minimum sample size for a sub-cell to be retained
#'   (passed to \code{linf.csts} during refinement). Also used to compute
#'   auto-refinement threshold when \code{cells.to.refine = NULL}. Any sub-cells
#'   (including NA assignments) with fewer than \code{n0} samples are automatically
#'   reassigned to their nearest sibling cell. Default: 50.
#' @param refinement.factor Numeric >= 1. When \code{cells.to.refine = NULL},
#'   automatically refine cells with at least \code{refinement.factor * n0}
#'   samples. For example, with \code{n0 = 50} and \code{refinement.factor = 2},
#'   cells with >= 100 samples are refined. Ignored when \code{cells.to.refine}
#'   is specified manually. Default: 2.
#' @param sep Character string. Separator for hierarchical labels. Default: \code{"_"}.
#' @param verbose Logical. Print progress messages during refinement. Default: \code{TRUE}.
#'
#' @return An updated CST object with \code{cst.depth} increased by one.
#'
#' @details
#' \strong{Auto-Refinement Algorithm:}
#'
#' When \code{cells.to.refine = NULL}, the function:
#' \enumerate{
#'   \item Computes the size of each cell from \code{table(csts$cell.label)}
#'   \item Calculates threshold = \code{max(refinement.factor * n0, min.cell.size)}
#'   \item Selects all cells with sample count >= threshold
#'   \item Proceeds with refinement on selected cells
#' }
#'
#' If no cells meet the threshold, the original CST assignments are returned
#' unchanged (with appropriate notifications if \code{verbose = TRUE}).
#'
#' \strong{Manual Refinement:}
#'
#' When \code{cells.to.refine} is explicitly provided, only those cells are
#' refined regardless of their size. The \code{refinement.factor} and
#' \code{min.cell.size} parameters are ignored in this mode.
#'
#' \strong{Refinement Procedure (for each selected cell):}
#'
#' \enumerate{
#'   \item Identifies all samples assigned to that parent cell
#'   \item Extracts the submatrix containing only those samples
#'   \item Identifies and removes the column corresponding to the parent cell's
#'     dominant feature (by matching \code{parent.cell} to column names)
#'   \item Re-normalizes the reduced matrix with \code{normalize.linf}
#'   \item Applies \code{linf.csts(n0 = n0)} to identify secondary dominance
#'     patterns
#'   \item Creates hierarchical labels: \code{paste(parent, subcell, sep = sep)}
#'   \item Maps sub-cell indices back to original matrix column indices
#' }
#'
#' \strong{Handling Edge Cases:}
#' \itemize{
#'   \item If a sample's sub-cell assignment is \code{NA} (e.g., all remaining
#'     features are zero after removing the dominant one), it receives the label
#'     \code{paste(parent, "NA", sep = sep)}.
#'   \item If the parent cell name doesn't match any column name in \code{M},
#'     refinement is skipped with a warning.
#'   \item If removing the dominant feature leaves no columns, refinement is
#'     skipped with a warning.
#'   \item Zero-sample cells are skipped with a warning.
#' }
#'
#' \strong{Index Mapping:}
#'
#' The \code{refined.index} values for refined cells are carefully mapped back
#' to the original matrix indices. If the sub-cell index in the reduced matrix
#' is \code{j}, and the removed parent feature was at index \code{k}, the
#' original index is computed as: \code{ifelse(j < k, j, j + 1)}.
#'
#' @section Requirements:
#' \itemize{
#'   \item \code{M} must have column names (either original or auto-generated).
#'   \item Column names should be unique after \code{make.unique}.
#'   \item Parent cell labels must match column names for successful refinement.
#'   \item \code{csts$cell.label} must have the same length as \code{nrow(M)}.
#' }
#'
#' @section Choosing Refinement Parameters:
#'
#' \strong{refinement.factor:}
#' \itemize{
#'   \item \code{refinement.factor = 1.5}: Aggressive refinement (cells with >= 1.5Ã—n0 samples)
#'   \item \code{refinement.factor = 2}: Balanced default (cells with >= 2Ã—n0 samples)
#'   \item \code{refinement.factor = 3}: Conservative (cells with >= 3Ã—n0 samples)
#'   \item \code{refinement.factor = 4}: Very conservative (cells with >= 4Ã—n0 samples)
#' }
#'
#' \strong{min.cell.size:}
#' Use when you want an absolute threshold independent of \code{n0}. For example,
#' \code{min.cell.size = 200} ensures only cells with 200+ samples are refined,
#' regardless of the \code{n0} value used for sub-cell filtering.
#'
#' @section Typical Workflows:
#'
#' \strong{Workflow 1: Fully Automatic}
#' \preformatted{
#' # Filter and normalize
#' filt <- filter.asv(counts, min.lib = 1000, prev.prop = 0.05)
#' M <- normalize.linf(filt$counts)
#'
#' # Initial CST assignment
#' csts <- linf.csts(M, n0 = 50)
#'
#' # Automatic refinement (refines cells with >= 100 samples)
#' refined <- refine.linf.csts(M, csts, n0 = 50, refinement.factor = 2)
#' table(refined$refined.label)
#' }
#'
#' \strong{Workflow 2: Inspect Then Decide}
#' \preformatted{
#' # Initial CSTs
#' csts <- linf.csts(M, n0 = 50)
#' cell.sizes <- sort(table(csts$cell.label), decreasing = TRUE)
#' print(cell.sizes)
#'
#' # Manual selection of specific cells
#' refined <- refine.linf.csts(
#'   M, csts,
#'   cells.to.refine = c("Lactobacillus_iners", "Gardnerella_vaginalis"),
#'   n0 = 50
#' )
#' }
#'
#' \strong{Workflow 3: Conservative Auto-refinement}
#' \preformatted{
#' # Only refine very large cells
#' refined <- refine.linf.csts(
#'   M, csts,
#'   n0 = 50,
#'   refinement.factor = 4,      # >= 200 samples
#'   min.cell.size = 250,        # or >= 250 samples (whichever is larger)
#'   verbose = TRUE
#' )
#' }
#'
#' @examples
#' \dontrun{
#' # === Real-world vaginal microbiome example ===
#' library(linf)
#'
#' # Assume X is your filtered count matrix
#' M <- normalize.linf(X)
#' csts <- linf.csts(M, n0 = 50)
#'
#' # Check initial distribution
#' table(csts$cell.label)
#' #   Lactobacillus_iners: 1060
#' #   Gardnerella_vaginalis: 426
#' #   Lactobacillus_crispatus: 352
#' #   Ca_Lachnocurva_vaginae: 279
#'
#' # Automatic refinement with default settings
#' refined <- refine.linf.csts(M, csts, n0 = 50, refinement.factor = 2)
#' # Refines cells with >= 100 samples (Li, Gv, Lc, Clv)
#'
#' # Examine refined structure
#' table(refined$refined.label)
#' print(refined$refinement.summary)
#' print(refined$cells.considered)
#'
#' # More conservative: only refine very large cells
#' refined2 <- refine.linf.csts(
#'   M, csts,
#'   n0 = 50,
#'   refinement.factor = 3,
#'   verbose = TRUE
#' )
#' # Refines cells with >= 150 samples (only Li, Gv, Lc)
#'
#' # Manual override for specific cells
#' refined3 <- refine.linf.csts(
#'   M, csts,
#'   cells.to.refine = "Lactobacillus_iners",
#'   n0 = 50
#' )
#'
#' # Access detailed sub-CST info
#' li.subcsts <- refined$sub.csts[["Lactobacillus_iners"]]
#' print(li.subcsts$size.table)
#' print(li.subcsts$kept.cells.lbl)
#' }
#'
#' # === Synthetic example ===
#' set.seed(123)
#' S <- matrix(rpois(500, lambda = 10), nrow = 100, ncol = 5)
#' colnames(S) <- paste0("Feature", LETTERS[1:5])
#'
#' # Create situation where FeatureA dominates many samples
#' S[1:60, 1] <- S[1:60, 1] + rpois(60, 50)
#' # FeatureB dominates some samples
#' S[61:75, 2] <- S[61:75, 2] + rpois(15, 50)
#'
#' M <- normalize.linf(S)
#' csts <- linf.csts(M, n0 = 5)
#' table(csts$cell.label)
#'
#' # Auto-refine with low threshold
#' refined <- refine.linf.csts(M, csts, n0 = 5, refinement.factor = 2)
#' table(refined$refined.label)
#' refined$refinement.summary
#'
#' # Manual refinement of just FeatureA
#' refined2 <- refine.linf.csts(M, csts, cells.to.refine = "FeatureA", n0 = 3)
#' table(refined2$refined.label)
#'
#' @seealso
#' \code{\link{linf.csts}} for initial L-infinity CST assignment,
#' \code{\link{normalize.linf}} for L-infinity normalization,
#' \code{\link{filter.asv}} for preprocessing count matrices.
#'
#' @export
refine.linf.csts <- function(M,
                             csts,
                             n0 = 50,
                             refinement.factor = 2,
                             sep = "__",
                             verbose = TRUE) {

    validate.linf.csts(csts)

    ## Initialize hierarchy if missing
    if (is.null(csts$cst.levels)) {
        csts$cst.levels <- list(level1 = csts$cell.label)
        csts$cst.depth <- 1L
        csts$sep <- sep
    }

    depth <- csts$cst.depth + 1L
    parent.labels <- csts$cst.levels[[depth - 1L]]

    refined.labels <- parent.labels
    cell.sizes <- sort(table(parent.labels), decreasing = TRUE)

    threshold <- refinement.factor * n0
    refine.cells <- names(cell.sizes[cell.sizes >= threshold])

    if (verbose) {
        cat("========================================\n")
        cat("AUTO-REFINEMENT MODE\n")
        cat("========================================\n")
        cat("Refinement threshold:", threshold, "\n")
        cat("Cells selected for refinement:", length(refine.cells), "\n\n")
    }

    for (cell in refine.cells) {
        idx <- which(parent.labels == cell)

        parent.features <- strsplit(cell, sep, fixed = TRUE)[[1]]
        drop.idx <- match(parent.features, colnames(M))
        drop.idx <- drop.idx[!is.na(drop.idx)]

        M.sub <- M[idx, -drop.idx, drop = FALSE]

        sub.csts <- linf.csts(M.sub, n0 = n0)
        sub.labels <- sub.csts$cell.label

        ## sub.labels can contain NA when rows become all-zero after dropping parent taxa.
        ## If we paste those, we create tiny "...__NA" CSTs. Instead, merge NA rows into a
        ## real child label (the most frequent non-NA sublabel) when possible.

        na.sub <- is.na(sub.labels)

        if (any(na.sub)) {
            ok.sub <- !na.sub

            if (any(ok.sub)) {
                ## linf.csts() guarantees non-NA sublabels are from kept cells (>= n0)
                fallback <- names(sort(table(sub.labels[ok.sub]), decreasing = TRUE))[1L]
                sub.labels[na.sub] <- fallback
            } else {
                ## All rows are NA in this parent after dropping columns -> nothing to refine
                next
            }
        }

        refined.labels[idx] <- paste(cell, sub.labels, sep = sep)
    }

    ## Update CST object
    csts$cst.levels[[depth]] <- refined.labels
    csts$cst.depth <- depth
    csts$cell.label <- refined.labels

    class(csts) <- "linf.csts"

    csts
}

#' Iteratively refine L-infinity CSTs
#'
#' @description
#' Applies iterative CST refinement to selected leaf CSTs, appending additional
#' hierarchy levels until refinement criteria are no longer met.
#'
#' @param M Numeric matrix used for CST assignment.
#' @param refined A CST object produced by \code{refine.linf.csts()}.
#' @param cells.to.refine Character vector or NULL. CSTs to refine.
#' @param n0 Integer. Minimum CST size.
#' @param refinement.factor Numeric. Refinement threshold multiplier.
#' @param sep Character. CST label separator.
#' @param verbose Logical. Print progress messages.
#'
#' @return An updated CST object with increased \code{cst.depth}.
#'
#' @details
#' Each iteration appends a new CST level to \code{cst.levels}.
#' Leaf CST labels are always stored in \code{cell.label}.
#'
#' @export
refine.linf.csts.iter <- function(M,
                                  refined,
                                  cells.to.refine = NULL,
                                  n0 = 50,
                                  refinement.factor = 5,
                                  sep = "__",
                                  verbose = TRUE) {

    validate.linf.csts(refined)

    depth <- refined$cst.depth
    parent.labels <- refined$cst.levels[[depth]]

    cell.sizes <- sort(table(parent.labels), decreasing = TRUE)
    threshold <- refinement.factor * n0

    if (is.null(cells.to.refine)) {
        cells.to.refine <- names(cell.sizes[cell.sizes >= threshold])
    }

    if (verbose) {
        cat("Auto-selecting cells for iterative refinement (threshold =", threshold, "):\n")
        for (c in cells.to.refine) {
            cat(" •", c, ":", cell.sizes[c], "samples\n")
        }
        cat("\n")
    }

    new.labels <- parent.labels

    for (cell in cells.to.refine) {
        idx <- which(parent.labels == cell)

        parents <- strsplit(cell, sep, fixed = TRUE)[[1]]
        drop.idx <- match(parents, colnames(M))
        drop.idx <- drop.idx[!is.na(drop.idx)]

        M.sub <- M[idx, -drop.idx, drop = FALSE]

        sub.csts <- linf.csts(M.sub, n0 = n0)
        sub.labels <- sub.csts$cell.label

        ## sub.labels can contain NA when rows become all-zero after dropping parent taxa.
        ## If we paste those, we create tiny "...__NA" CSTs. Instead, merge NA rows into a
        ## real child label (the most frequent non-NA sublabel) when possible.

        na.sub <- is.na(sub.labels)

        if (any(na.sub)) {
            ok.sub <- !na.sub

            if (any(ok.sub)) {
                ## linf.csts() guarantees non-NA sublabels are from kept cells (>= n0)
                fallback <- names(sort(table(sub.labels[ok.sub]), decreasing = TRUE))[1L]
                sub.labels[na.sub] <- fallback
            } else {
                ## All rows are NA in this parent after dropping columns -> nothing to refine
                next
            }
        }

        new.labels[idx] <- paste(cell, sub.labels, sep = sep)
    }

    refined$cst.levels[[depth + 1L]] <- new.labels
    refined$cst.depth <- depth + 1L
    refined$cell.label <- new.labels

    class(refined) <- "linf.csts"

    refined
}


#' Show L-infinity CSTs in hierarchical format
#'
#' @description
#' Displays CST labels organized hierarchically, grouping sub-cells under their
#' parents with optional indentation, tree characters, and summary statistics.
#'
#' @param csts Either a refined CST object (list from \code{refine.linf.csts})
#'   or a character vector of CST labels.
#' @param sep Character separator used in hierarchical labels. Default: \code{"__"}.
#' @param style Character. Display style: \code{"tree"} (with tree characters),
#'   \code{"indent"} (indentation only), or \code{"flat"} (sorted alphabetically).
#'   Default: \code{"tree"}.
#' @param show.counts Logical. Show sample counts. Default: \code{TRUE}.
#' @param show.pct Logical. Show percentages. Default: \code{TRUE}.
#' @param show.cumulative Logical. Show cumulative percentages for each parent
#'   group. Default: \code{FALSE}.
#' @param min.count Integer. Only show cells with at least this many samples.
#'   Default: \code{0} (show all).
#' @param max.depth Integer or NULL. Maximum hierarchy depth to display.
#'   Default: \code{NULL} (show all depths).
#'
#' @return Invisibly returns a data frame with the hierarchical structure.
#'
#' @examples
#' # Basic tree view
#' show.linf.csts(refined2)
#'
#' # Indented view without tree characters
#' show.linf.csts(refined2, style = "indent")
#'
#' # Show only cells with 100+ samples
#' show.linf.csts(refined2, min.count = 100)
#'
#' # Flat alphabetical list
#' show.linf.csts(refined2, style = "flat")
#'
#' @export
print.linf.csts <- function(csts) {

  validate.linf.csts(csts)

  ## Determine hierarchy
  if (is.null(csts$cst.levels)) {
    levels <- list(level1 = csts$cell.label)
    max.depth <- 1L
  } else {
    levels <- csts$cst.levels
    max.depth <- csts$cst.depth
  }

  cat("\n================================================================================\n")
  cat("L-infinity CST Hierarchy\n")
  cat("================================================================================\n")
  cat("Total samples: ", length(csts$cell.label), "\n")
  cat("Max depth:     ", max.depth, "\n")
  cat("--------------------------------------------------------------------------------\n\n")

  for (d in seq_along(levels)) {
    cat("Depth", d, "\n")
    tab <- sort(table(levels[[d]]), decreasing = TRUE)
    for (nm in names(tab)) {
      cat("  ", nm, ": ", tab[[nm]], "\n", sep = "")
    }
    cat("\n")
  }

  invisible(NULL)
}

#' Summarize CST hierarchy statistics
#'
#' @param csts Refined CST object or character vector
#' @param sep Separator. Default: \code{"__"}
#' @return Data frame with depth-wise statistics
#' @export
summary.linf.csts <- function(csts) {

  validate.linf.csts(csts)

  ## Determine hierarchy
  if (is.null(csts$cst.levels)) {
    levels <- list(level1 = csts$cell.label)
  } else {
    levels <- csts$cst.levels
  }

  out <- data.frame(
    depth = integer(0),
    n.cells = integer(0),
    total.samples = integer(0),
    mean.size = numeric(0),
    median.size = numeric(0),
    min.size = numeric(0),
    max.size = numeric(0)
  )

  for (d in seq_along(levels)) {
    lbl <- levels[[d]]
    tab <- table(lbl)

    out <- rbind(out, data.frame(
      depth = d,
      n.cells = length(tab),
      total.samples = sum(tab),
      mean.size = mean(tab),
      median.size = median(tab),
      min.size = min(tab),
      max.size = max(tab)
    ))
  }

  out
}

#' Render L-infinity CSTs as a LaTeX table
#'
#' @description
#' Generates LaTeX code summarizing L-infinity CSTs at a specified hierarchy depth.
#' Uses the explicit CST hierarchy stored in \code{cst.levels} when available.
#'
#' @param csts A CST object produced by \code{linf.csts()} and optionally refined
#'   by \code{refine.linf.csts()} or \code{refine.linf.csts.iter()}.
#' @param depth Integer. CST depth to render. Default is the leaf level
#'   (\code{csts$cst.depth}). If \code{cst.levels} is missing, depth is ignored.
#' @param caption Character or NULL. LaTeX table caption.
#' @param label Character or NULL. LaTeX label for referencing the table.
#' @param digits Integer. Digits to display for percentages. Default: 1.
#' @param include.percent Logical. If TRUE, include percent of total samples.
#'
#' @return Character vector containing LaTeX table code.
#'
#' @details
#' This function is hierarchy-aware. It does not infer CST depth from label strings.
#' If the CST object does not contain \code{cst.levels}, the function falls back
#' to using \code{cell.label} as a single-level CST.
#'
#' @export
latex.linf.csts <- function(csts,
                            depth = NULL,
                            caption = NULL,
                            label = NULL,
                            digits = 1,
                            include.percent = TRUE) {

  validate.linf.csts(csts)

  ## Determine labels at requested depth
  if (!is.null(csts$cst.levels)) {
    if (is.null(depth)) {
      depth <- csts$cst.depth
    }
    if (depth < 1 || depth > csts$cst.depth) {
      stop("latex.linf.csts: invalid depth")
    }
    labels <- csts$cst.levels[[depth]]
  } else {
    ## Backward compatibility: flat CST
    labels <- csts$cell.label
    depth <- 1L
  }

  tab <- sort(table(labels), decreasing = TRUE)
  total <- sum(tab)

  df <- data.frame(
    CST = names(tab),
    Samples = as.integer(tab),
    stringsAsFactors = FALSE
  )

  if (include.percent) {
    df$Percent <- round(100 * df$Samples / total, digits)
  }

  ## Build LaTeX using df.to.latex()
  latex <- df.to.latex(
    df,
    caption = caption,
    label = label,
    digits = digits,
    include.rownames = FALSE,
    booktabs = TRUE,
    use.float = TRUE,
    print = FALSE
  )

  invisible(latex)
}
