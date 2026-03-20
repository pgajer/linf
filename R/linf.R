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

  ## Identify rows with meaningful L-infinity norm
  keep <- m > tol

  if (any(keep)) {
    X[keep, ] <- X[keep, , drop = FALSE] / m[keep]
  }

  X
}

#' L-infinity cell assignment
#'
#' @description
#' Assigns each row to the column achieving its maximum (L-infinity cell).
#'
#' For each sample (row) of a nonnegative matrix, identifies the \eqn{L^\infty}
#' cell as the column with the maximum value. Ties are broken by the first
#' maximum (as in \code{max.col(..., ties.method = "first")}). Rows that are all
#' zero are assigned \code{NA}.
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
#' @param tie.method Character. How to resolve ties when assigning L-infinity cells.
#' @param return.value Logical. If `TRUE`, include a `value` vector with row maxima.
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
                       feature.ids = NULL,
                       feature.labels = NULL,
                       tie.method = c("first", "random", "error"),
                       return.value = FALSE) {

  tie.method <- match.arg(tie.method)

  X <- as.matrix(S)
  storage.mode(X) <- "numeric"

  if (any(!is.finite(X))) stop("linf.cells: non-finite entries found")
  if (any(X < 0, na.rm = TRUE)) stop("linf.cells: negative entries found")
  if (nrow(X) == 0L || ncol(X) == 0L)
    stop("linf.cells: matrix has zero rows or columns")

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
      stop("linf.cells: tie encountered and tie.method = 'error'")
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

#' Truncated L-infinity CSTs with configurable low-frequency handling
#'
#' @description
#' Computes L-infinity cells (dominant feature per sample) and then applies a
#' minimum size threshold \code{n0}. Cells with fewer than \code{n0} samples are
#' handled according to \code{low.freq.policy}:
#' \itemize{
#'   \item \code{"rare"}: collapse all low-frequency cells into \code{rare.label}.
#'   \item \code{"absorb"}: reassign each low-frequency sample to the kept cell
#'     with the largest value among the kept set (ties handled by \code{tie.method}).
#' }
#'
#' @param S Numeric matrix (samples x features), typically L-infinity relatives.
#' @param feature.ids Optional character vector of stable feature identifiers,
#'   length \code{ncol(S)}.
#' @param feature.labels Optional character vector of display labels, length
#'   \code{ncol(S)}.
#' @param n0 Integer >= 1. Minimum size for a cell to be kept.
#' @param low.freq.policy Character. One of \code{"rare"} or \code{"absorb"}.
#'   Default: \code{"rare"}.
#' @param rare.label Character scalar used when \code{low.freq.policy = "rare"}.
#'   Default: \code{"RARE_DOMINANT"}.
#' @param tie.method Character. Tie handling passed to \code{linf.cells()} and used
#'   during absorb reassignment ("first", "random", "error").
#' @param return.diagnostics Logical. If TRUE, return reassignment diagnostics.
#' @param return.landmarks Logical. If TRUE, attach a depth-1 landmark summary
#'   computed by \code{\link{linf.landmarks}}.
#' @param landmark.types Character vector of landmark types passed to
#'   \code{\link{linf.landmarks}} when \code{return.landmarks = TRUE}.
#' @param landmark.view Character. Landmark view passed to
#'   \code{\link{linf.landmarks}} when \code{return.landmarks = TRUE}.
#'
#' @return List with:
#'   \itemize{
#'     \item \code{cell.index}, \code{cell.id}, \code{cell.label}: active labeling per \code{low.freq.policy}
#'     \item \code{cell.index.rare}, \code{cell.id.rare}, \code{cell.label.rare}
#'     \item \code{cell.index.absorb}, \code{cell.id.absorb}, \code{cell.label.absorb}
#'     \item \code{kept.cells.idx}, \code{kept.cells.id}, \code{kept.cells.lbl}
#'     \item \code{raw.index}, \code{raw.id}, \code{raw.label}
#'     \item \code{feature.ids}, \code{feature.labels}
#'     \item \code{size.table}, \code{size.table.id}
#'     \item \code{n0}, \code{low.freq.policy}, \code{rare.label}
#'     \item \code{diagnostics} (if \code{return.diagnostics = TRUE})
#'     \item \code{landmarks} (if \code{return.landmarks = TRUE})
#'   }
#'
#' @export
linf.csts <- function(S,
                      feature.ids = NULL,
                      feature.labels = NULL,
                      n0 = 50,
                      low.freq.policy = c("rare", "absorb"),
                      rare.label = "RARE_DOMINANT",
                      tie.method = c("first", "random", "error"),
                      return.diagnostics = FALSE,
                      return.landmarks = FALSE,
                      landmark.types = c("endpoint.max", "endpoint.min"),
                      landmark.view = c("active", "rare", "absorb")) {

    low.freq.policy <- match.arg(low.freq.policy)
    tie.method <- match.arg(tie.method)
    landmark.view <- match.arg(landmark.view)

    if (!is.numeric(n0) || length(n0) != 1L || n0 < 1 || n0 %% 1 != 0) {
        stop("linf.csts: n0 must be integer >= 1")
    }
    if (!is.character(rare.label) || length(rare.label) != 1L || !nzchar(rare.label)) {
        stop("linf.csts: rare.label must be a non-empty character scalar")
    }

    X <- as.matrix(S)
    storage.mode(X) <- "numeric"

    if (any(!is.finite(X))) stop("linf.csts: non-finite entries found")
    if (any(X < 0, na.rm = TRUE)) stop("linf.csts: negative entries found")
    if (nrow(X) == 0L || ncol(X) == 0L)
        stop("linf.csts: matrix has zero rows or columns")

    meta <- resolve.linf.feature.meta(X, feature.ids = feature.ids, feature.labels = feature.labels)
    fid <- meta$feature.ids
    lev <- meta$feature.labels

    raw <- linf.cells(X,
                      feature.ids = fid,
                      feature.labels = lev,
                      tie.method = tie.method)
    tab <- sort(table(raw$label[!is.na(raw$label)]), decreasing = TRUE)
    tab.id <- sort(table(raw$id[!is.na(raw$id)]), decreasing = TRUE)

    kept.lbl <- names(tab[tab >= n0])
    kept.idx <- match(kept.lbl, lev)
    kept.id <- fid[kept.idx]

    n <- nrow(X)

    ## Rare-policy labels: keep only cells with size >= n0; everything else -> rare.label
    is.kept <- !is.na(raw$label) & (raw$label %in% kept.lbl)

    cell.idx.rare <- raw$index
    cell.id.rare <- raw$id
    cell.lbl.rare <- raw$label
    cell.idx.rare[!is.kept] <- NA_integer_
    cell.id.rare[!is.kept] <- rare.label
    cell.lbl.rare[!is.kept] <- rare.label

    ## Absorb-policy labels: reassign low-frequency (and zero-row) samples into kept cells
    cell.idx.absorb <- raw$index
    cell.id.absorb <- raw$id
    cell.lbl.absorb <- raw$label

    reassigned <- logical(n)
    reassigned.from <- rep(NA_character_, n)
    reassigned.to   <- rep(NA_character_, n)

    if (length(kept.idx) > 0L) {

        ## Fallback target for degenerate cases (e.g., all kept values are zero)
        fallback.idx <- kept.idx[1L]
        fallback.lbl <- lev[fallback.idx]

        ## Absorb: (i) low-frequency raw labels, (ii) raw NA labels (e.g., all-zero rows)
        to.absorb <- which(!is.kept)

        if (length(to.absorb)) {
            new.idx <- apply(X[to.absorb, , drop = FALSE], 1, function(r) {
                kvals <- r[kept.idx]
                m <- max(kvals)

                ## If there is no positive evidence among kept taxa, avoid arbitrary ties
                if (!is.finite(m) || m <= 0) return(fallback.idx)

                j <- kept.idx[kvals == m]
                if (length(j) == 1L) return(j)
                if (tie.method == "first") return(j[1L])
                if (tie.method == "random") return(sample(j, 1L))
                stop("linf.csts: tie during reassignment and tie.method = 'error'")
            })

            reassigned[to.absorb] <- TRUE
            reassigned.from[to.absorb] <- raw$label[to.absorb]
            reassigned.to[to.absorb]   <- lev[new.idx]

            cell.idx.absorb[to.absorb] <- new.idx
            cell.id.absorb[to.absorb] <- fid[new.idx]
            cell.lbl.absorb[to.absorb] <- lev[new.idx]
        }

    } else {

        ## No kept cells at this n0:
        ## - rare-policy: everyone is rare.label (already set above)
        ## - absorb-policy: undefined; keep NA labels
        cell.idx.absorb[] <- NA_integer_
        cell.id.absorb[] <- NA_character_
        cell.lbl.absorb[] <- NA_character_
    }

    ## Select active labeling
    if (low.freq.policy == "rare") {
        cell.idx <- cell.idx.rare
        cell.id <- cell.id.rare
        cell.lbl <- cell.lbl.rare
    } else {
        cell.idx <- cell.idx.absorb
        cell.id <- cell.id.absorb
        cell.lbl <- cell.lbl.absorb
    }

    names(cell.idx) <- rownames(X)
    names(cell.id) <- rownames(X)
    names(cell.lbl) <- rownames(X)

    names(cell.idx.rare) <- rownames(X)
    names(cell.id.rare) <- rownames(X)
    names(cell.lbl.rare) <- rownames(X)

    names(cell.idx.absorb) <- rownames(X)
    names(cell.id.absorb) <- rownames(X)
    names(cell.lbl.absorb) <- rownames(X)

    out <- list(
        cell.index       = cell.idx,
        cell.id          = cell.id,
        cell.label       = cell.lbl,
        cell.index.rare  = cell.idx.rare,
        cell.id.rare     = cell.id.rare,
        cell.label.rare  = cell.lbl.rare,
        cell.index.absorb = cell.idx.absorb,
        cell.id.absorb   = cell.id.absorb,
        cell.label.absorb = cell.lbl.absorb,
        kept.cells.idx   = kept.idx,
        kept.cells.id    = kept.id,
        kept.cells.lbl   = kept.lbl,
        raw.index        = raw$index,
        raw.id           = raw$id,
        raw.label        = raw$label,
        size.table       = tab,
        size.table.id    = tab.id,
        feature.ids      = fid,
        feature.labels   = lev,
        n0               = as.integer(n0),
        low.freq.policy  = low.freq.policy,
        rare.label       = rare.label
    )

    if (return.diagnostics) {
        out$diagnostics <- list(
            reassigned      = reassigned,
            reassigned.from = reassigned.from,
            reassigned.to   = reassigned.to
        )
    }

    out$cst.id.levels <- list(level1 = out$cell.id)
    out$cst.id.levels.rare <- list(level1 = out$cell.id.rare)
    out$cst.id.levels.absorb <- list(level1 = out$cell.id.absorb)
    out$cst.levels <- list(level1 = out$cell.label)
    out$cst.levels.rare <- list(level1 = out$cell.label.rare)
    out$cst.levels.absorb <- list(level1 = out$cell.label.absorb)
    out$cst.depth <- 1L

    class(out) <- "linf.csts"

    if (isTRUE(return.landmarks)) {
        out$landmarks <- linf.landmarks(
            X,
            out,
            depth = 1L,
            view = landmark.view,
            landmark.types = landmark.types,
            tie.method = tie.method
        )
    }

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

#' Refine L-infinity CST hierarchy by one level
#'
#' @description
#' Selects large leaf cells and refines them by dropping the dominant feature(s)
#' encoded in the parent label path and re-applying \code{\link{linf.csts}} to
#' the remaining features. The resulting child labels are appended to the parent
#' label using \code{sep}.
#'
#' Low-frequency child cells are handled by \code{low.freq.policy}. When
#' \code{low.freq.policy = "rare"}, rare buckets at depth >= 2 become
#' parent-prefixed automatically via the hierarchical \code{paste(parent, child, sep = sep)}.
#'
#' @param M Numeric matrix (samples x features) used for refinement. Column names
#'   should match feature labels used in CST names.
#' @param csts A \code{"linf.csts"} object.
#' @param n0 Integer >= 1. Minimum size for a child cell to be kept (passed to \code{linf.csts}).
#' @param refinement.factor Numeric > 0. Auto-refine parent cells with size >= \code{refinement.factor * n0}.
#' @param sep Character scalar used to concatenate hierarchical labels.
#' @param low.freq.policy Character. One of \code{"rare"} or \code{"absorb"}. Default: \code{"rare"}.
#' @param rare.label Character scalar for rare buckets when \code{low.freq.policy = "rare"}.
#' @param verbose Logical. If TRUE, print progress information.
#'
#' @return Updated \code{"linf.csts"} object with \code{cst.depth} increased by one and
#'   updated \code{cell.label}. Policy-specific views are stored in
#'   \code{cell.label.rare} and \code{cell.label.absorb}.
#'
#' @export
refine.linf.csts <- function(M,
                             csts,
                             n0 = 50,
                             refinement.factor = 2,
                             sep = "__",
                             low.freq.policy = c("rare", "absorb"),
                             rare.label = "RARE_DOMINANT",
                             verbose = TRUE) {

    validate.linf.csts(csts)

    low.freq.policy <- match.arg(low.freq.policy)

    if (!is.numeric(n0) || length(n0) != 1L || n0 < 1 || n0 %% 1 != 0) {
        stop("refine.linf.csts: n0 must be integer >= 1")
    }
    if (!is.numeric(refinement.factor) || length(refinement.factor) != 1L ||
        !is.finite(refinement.factor) || refinement.factor <= 0) {
        stop("refine.linf.csts: refinement.factor must be a finite numeric > 0")
    }
    if (!is.character(sep) || length(sep) != 1L || !nzchar(sep)) {
        stop("refine.linf.csts: sep must be a non-empty character scalar")
    }
    if (!is.character(rare.label) || length(rare.label) != 1L || !nzchar(rare.label)) {
        stop("refine.linf.csts: rare.label must be a non-empty character scalar")
    }

    ## Initialize hierarchy if missing
    if (is.null(csts$cst.levels)) {
        csts$cst.levels <- list(level1 = csts$cell.label)
        csts$cst.depth <- 1L
        csts$sep <- sep

        ## Ensure both policy-specific views exist at depth 1
        if (is.null(csts$cell.label.rare)) {
            csts$cell.label.rare <- csts$cell.label
        }
        if (is.null(csts$cell.label.absorb)) {
            csts$cell.label.absorb <- csts$cell.label
        }

        csts$cst.levels.rare <- list(level1 = csts$cell.label.rare)
        csts$cst.levels.absorb <- list(level1 = csts$cell.label.absorb)
    }

    if (is.null(csts$cst.id.levels)) {
        csts$cst.id.levels <- list(level1 = csts$cell.id %||% csts$cell.label)
    }

    if (is.null(csts$feature.ids)) {
        csts$feature.ids <- colnames(M)
    }
    if (is.null(csts$feature.labels)) {
        csts$feature.labels <- csts$feature.ids
    }

    ## Backward-compatible defaults if stored views are missing
    if (is.null(csts$cst.levels.rare)) csts$cst.levels.rare <- csts$cst.levels
    if (is.null(csts$cst.levels.absorb)) csts$cst.levels.absorb <- csts$cst.levels
    if (is.null(csts$cst.id.levels.rare)) csts$cst.id.levels.rare <- csts$cst.id.levels
    if (is.null(csts$cst.id.levels.absorb)) csts$cst.id.levels.absorb <- csts$cst.id.levels

    depth <- csts$cst.depth + 1L
    parent.ids <- csts$cst.id.levels[[depth - 1L]]
    parent.ids.rare <- csts$cst.id.levels.rare[[depth - 1L]]
    parent.ids.absorb <- csts$cst.id.levels.absorb[[depth - 1L]]
    parent.labels <- csts$cst.levels[[depth - 1L]]

    parent.labels.rare <- csts$cst.levels.rare[[depth - 1L]]
    parent.labels.absorb <- csts$cst.levels.absorb[[depth - 1L]]

    refined.ids.rare <- parent.ids.rare
    refined.ids.absorb <- parent.ids.absorb
    refined.labels.rare <- parent.labels.rare
    refined.labels.absorb <- parent.labels.absorb

    cell.sizes <- sort(table(parent.ids), decreasing = TRUE)

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
        idx <- which(parent.ids == cell)
        parent.id.rare <- parent.ids.rare[idx[1L]]
        parent.id.absorb <- parent.ids.absorb[idx[1L]]
        parent.label.rare <- parent.labels.rare[idx[1L]]
        parent.label.absorb <- parent.labels.absorb[idx[1L]]

        parent.features <- strsplit(cell, sep, fixed = TRUE)[[1]]
        drop.idx <- match(parent.features, csts$feature.ids)
        drop.idx <- drop.idx[!is.na(drop.idx)]

        ## If no parent features match columns (e.g., cell is a rare bucket), do not drop any columns.
        ## Note: x[, -integer(0)] selects *zero* columns, so we must handle this explicitly.
        if (length(drop.idx) == 0L) {
            M.sub <- M[idx, , drop = FALSE]
        } else if (length(drop.idx) >= ncol(M)) {
            ## Dropping all columns would yield an empty matrix; nothing to refine.
            next
        } else {
            M.sub <- M[idx, -drop.idx, drop = FALSE]
        }

        if (nrow(M.sub) == 0L || ncol(M.sub) == 0L) next

        sub.csts <- linf.csts(M.sub,
                             feature.ids = csts$feature.ids[-drop.idx],
                             feature.labels = csts$feature.labels[-drop.idx],
                             n0 = n0,
                             low.freq.policy = low.freq.policy,
                             rare.label = rare.label)

        sub.ids.rare <- sub.csts$cell.id.rare %||% sub.csts$cell.label.rare
        sub.ids.absorb <- sub.csts$cell.id.absorb %||% sub.csts$cell.label.absorb
        sub.labels.rare <- sub.csts$cell.label.rare
        sub.labels.absorb <- sub.csts$cell.label.absorb

        ## If refinement yields no kept sub-cells (all samples go to rare), skip
        if (all(sub.labels.rare == rare.label)) next

        refined.ids.rare[idx] <- paste(parent.id.rare, sub.ids.rare, sep = sep)
        refined.ids.absorb[idx] <- paste(parent.id.absorb, sub.ids.absorb, sep = sep)
        refined.labels.rare[idx] <- paste(parent.label.rare, sub.labels.rare, sep = sep)
        refined.labels.absorb[idx] <- paste(parent.label.absorb, sub.labels.absorb, sep = sep)
    }

    ## Update CST object (store both views; active view chosen by low.freq.policy)
    csts$cst.id.levels[[depth]] <- if (low.freq.policy == "rare") refined.ids.rare else refined.ids.absorb
    csts$cst.id.levels.rare[[depth]] <- refined.ids.rare
    csts$cst.id.levels.absorb[[depth]] <- refined.ids.absorb
    csts$cst.levels.rare[[depth]] <- refined.labels.rare
    csts$cst.levels.absorb[[depth]] <- refined.labels.absorb

    if (low.freq.policy == "rare") {
        refined.labels <- refined.labels.rare
    } else {
        refined.labels <- refined.labels.absorb
    }

    csts$cst.levels[[depth]] <- refined.labels
    csts$cst.depth <- depth
    csts$sep <- sep

    csts$cell.id.rare <- refined.ids.rare
    csts$cell.id.absorb <- refined.ids.absorb
    csts$cell.id <- csts$cst.id.levels[[depth]]
    csts$cell.label.rare <- refined.labels.rare
    csts$cell.label.absorb <- refined.labels.absorb
    csts$cell.label <- refined.labels

    csts$low.freq.policy <- low.freq.policy
    csts$rare.label <- rare.label
    csts$landmarks <- NULL

    class(csts) <- "linf.csts"

    csts
}

#' Iteratively refine L-infinity CSTs by one additional level
#'
#' @description
#' Appends one refinement level to an existing CST hierarchy produced by
#' \code{\link{refine.linf.csts}} or \code{\link{refine.linf.csts.iter}}.
#' Cells to refine can be provided explicitly or selected automatically based on
#' \code{refinement.factor * n0}.
#'
#' @param M Numeric matrix used for refinement.
#' @param refined A \code{"linf.csts"} object with an existing hierarchy.
#' @param cells.to.refine Character vector of leaf cells to refine, or NULL for auto-selection.
#' @param n0 Integer >= 1. Minimum size for a child cell to be kept.
#' @param refinement.factor Numeric > 0. Auto-selection threshold multiplier.
#' @param sep Character scalar used to concatenate hierarchical labels.
#' @param low.freq.policy Character. One of \code{"rare"} or \code{"absorb"}. Default: \code{"rare"}.
#' @param rare.label Character scalar for rare buckets when \code{low.freq.policy = "rare"}.
#' @param verbose Logical. If TRUE, print progress information.
#'
#' @return Updated \code{"linf.csts"} object with \code{cst.depth} increased by one.
#'
#' @export
refine.linf.csts.iter <- function(M,
                                  refined,
                                  cells.to.refine = NULL,
                                  n0 = 50,
                                  refinement.factor = 5,
                                  sep = "__",
                                  low.freq.policy = c("rare", "absorb"),
                                  rare.label = "RARE_DOMINANT",
                                  verbose = TRUE) {

    validate.linf.csts(refined)

    low.freq.policy <- match.arg(low.freq.policy)

    if (!is.numeric(n0) || length(n0) != 1L || n0 < 1 || n0 %% 1 != 0) {
        stop("refine.linf.csts.iter: n0 must be integer >= 1")
    }
    if (!is.numeric(refinement.factor) || length(refinement.factor) != 1L ||
        !is.finite(refinement.factor) || refinement.factor <= 0) {
        stop("refine.linf.csts.iter: refinement.factor must be a finite numeric > 0")
    }
    if (!is.character(sep) || length(sep) != 1L || !nzchar(sep)) {
        stop("refine.linf.csts.iter: sep must be a non-empty character scalar")
    }
    if (!is.character(rare.label) || length(rare.label) != 1L || !nzchar(rare.label)) {
        stop("refine.linf.csts.iter: rare.label must be a non-empty character scalar")
    }

    ## Ensure stored views exist (backward compatible)
    if (is.null(refined$cst.levels.rare)) refined$cst.levels.rare <- refined$cst.levels
    if (is.null(refined$cst.levels.absorb)) refined$cst.levels.absorb <- refined$cst.levels
    if (is.null(refined$cst.id.levels)) refined$cst.id.levels <- list(level1 = refined$cell.id %||% refined$cell.label)
    if (is.null(refined$cst.id.levels.rare)) refined$cst.id.levels.rare <- refined$cst.id.levels
    if (is.null(refined$cst.id.levels.absorb)) refined$cst.id.levels.absorb <- refined$cst.id.levels
    if (is.null(refined$feature.ids)) refined$feature.ids <- colnames(M)
    if (is.null(refined$feature.labels)) refined$feature.labels <- refined$feature.ids

    depth <- refined$cst.depth

    parent.ids <- refined$cst.id.levels[[depth]]
    parent.ids.rare <- refined$cst.id.levels.rare[[depth]]
    parent.ids.absorb <- refined$cst.id.levels.absorb[[depth]]
    parent.labels <- refined$cst.levels[[depth]]
    parent.labels.rare <- refined$cst.levels.rare[[depth]]
    parent.labels.absorb <- refined$cst.levels.absorb[[depth]]

    cell.sizes <- sort(table(parent.ids), decreasing = TRUE)
    threshold <- refinement.factor * n0

    if (is.null(cells.to.refine)) {
        cells.to.refine <- names(cell.sizes[cell.sizes >= threshold])
    }

    if (verbose) {
        cat("Auto-selecting cells for iterative refinement (threshold =", threshold, "):\n")
        for (c in cells.to.refine) {
            cat(" -", c, ":", cell.sizes[c], "samples\n")
        }
        cat("\n")
    }

    new.ids.rare <- parent.ids.rare
    new.ids.absorb <- parent.ids.absorb
    new.labels.rare <- parent.labels.rare
    new.labels.absorb <- parent.labels.absorb

    for (cell in cells.to.refine) {
        idx <- which(parent.ids == cell)
        parent.id.rare <- parent.ids.rare[idx[1L]]
        parent.id.absorb <- parent.ids.absorb[idx[1L]]
        parent.label.rare <- parent.labels.rare[idx[1L]]
        parent.label.absorb <- parent.labels.absorb[idx[1L]]

        parents <- strsplit(cell, sep, fixed = TRUE)[[1]]
        drop.idx <- match(parents, refined$feature.ids)
        drop.idx <- drop.idx[!is.na(drop.idx)]

        ## If no parent features match columns (e.g., cell is a rare bucket), do not drop any columns.
        ## Note: x[, -integer(0)] selects *zero* columns, so we must handle this explicitly.
        if (length(drop.idx) == 0L) {
            M.sub <- M[idx, , drop = FALSE]
        } else if (length(drop.idx) >= ncol(M)) {
            ## Dropping all columns would yield an empty matrix; nothing to refine.
            next
        } else {
            M.sub <- M[idx, -drop.idx, drop = FALSE]
        }

        if (nrow(M.sub) == 0L || ncol(M.sub) == 0L) next

        sub.csts <- linf.csts(M.sub,
                             feature.ids = refined$feature.ids[-drop.idx],
                             feature.labels = refined$feature.labels[-drop.idx],
                             n0 = n0,
                             low.freq.policy = low.freq.policy,
                             rare.label = rare.label)

        sub.ids.rare <- sub.csts$cell.id.rare %||% sub.csts$cell.label.rare
        sub.ids.absorb <- sub.csts$cell.id.absorb %||% sub.csts$cell.label.absorb
        sub.labels.rare <- sub.csts$cell.label.rare
        sub.labels.absorb <- sub.csts$cell.label.absorb

        ## If refinement yields no kept sub-cells (all samples go to rare), skip
        if (all(sub.labels.rare == rare.label)) next

        new.ids.rare[idx] <- paste(parent.id.rare, sub.ids.rare, sep = sep)
        new.ids.absorb[idx] <- paste(parent.id.absorb, sub.ids.absorb, sep = sep)
        new.labels.rare[idx] <- paste(parent.label.rare, sub.labels.rare, sep = sep)
        new.labels.absorb[idx] <- paste(parent.label.absorb, sub.labels.absorb, sep = sep)
    }

    refined$cst.id.levels[[depth + 1L]] <- if (low.freq.policy == "rare") new.ids.rare else new.ids.absorb
    refined$cst.id.levels.rare[[depth + 1L]] <- new.ids.rare
    refined$cst.id.levels.absorb[[depth + 1L]] <- new.ids.absorb
    refined$cst.levels.rare[[depth + 1L]] <- new.labels.rare
    refined$cst.levels.absorb[[depth + 1L]] <- new.labels.absorb

    if (low.freq.policy == "rare") {
        new.labels <- new.labels.rare
    } else {
        new.labels <- new.labels.absorb
    }

    refined$cst.levels[[depth + 1L]] <- new.labels
    refined$cst.depth <- depth + 1L

    refined$cell.id.rare <- new.ids.rare
    refined$cell.id.absorb <- new.ids.absorb
    refined$cell.id <- refined$cst.id.levels[[depth + 1L]]
    refined$cell.label.rare <- new.labels.rare
    refined$cell.label.absorb <- new.labels.absorb
    refined$cell.label <- new.labels

    refined$low.freq.policy <- low.freq.policy
    refined$rare.label <- rare.label
    refined$sep <- sep
    refined$landmarks <- NULL

    class(refined) <- "linf.csts"

    refined
}

#' Switch a CST hierarchy to the "absorb" view (collapse rare buckets)
#'
#' @description
#' Returns a copy of a \code{"linf.csts"} object with \code{cell.label} (and, if
#' present, \code{cst.levels}) replaced by the precomputed "absorb" labeling.
#' This does not recompute CSTs; it only switches between labelings already
#' stored in the object.
#'
#' @param csts A \code{"linf.csts"} object produced by \code{\link{linf.csts}}
#'   and optionally refined by \code{\link{refine.linf.csts}} /
#'   \code{\link{refine.linf.csts.iter}}.
#'
#' @return A \code{"linf.csts"} object using the "absorb" view.
#'
#' @export
collapse.rare <- function(csts) {

  validate.linf.csts(csts)

  if (is.null(csts$cell.label.absorb)) {
    stop("collapse.rare: no precomputed absorb labels found in object")
  }

  if (!is.null(csts$cell.index.absorb)) {
    csts$cell.index <- csts$cell.index.absorb
  }
  if (!is.null(csts$cell.id.absorb)) {
    csts$cell.id <- csts$cell.id.absorb
  }
  csts$cell.label <- csts$cell.label.absorb
  csts$low.freq.policy <- "absorb"

  if (!is.null(csts$cst.id.levels)) {
    if (is.null(csts$cst.id.levels.absorb)) {
      depth <- csts$cst.depth
      csts$cst.id.levels[[depth]] <- csts$cell.id %||% csts$cell.label
    } else {
      csts$cst.id.levels <- csts$cst.id.levels.absorb
    }
  }

  if (!is.null(csts$cst.levels)) {
    if (is.null(csts$cst.levels.absorb)) {
      ## Fall back to updating only the leaf level
      depth <- csts$cst.depth
      csts$cst.levels[[depth]] <- csts$cell.label.absorb
    } else {
      csts$cst.levels <- csts$cst.levels.absorb
    }
  }

  csts$landmarks <- NULL
  class(csts) <- "linf.csts"
  csts
}

#' Switch a CST hierarchy to the "rare" view (explicit rare buckets)
#'
#' @description
#' Returns a copy of a \code{"linf.csts"} object with \code{cell.label} (and, if
#' present, \code{cst.levels}) replaced by the precomputed "rare" labeling.
#' This does not recompute CSTs; it only switches between labelings already
#' stored in the object.
#'
#' @param csts A \code{"linf.csts"} object produced by \code{\link{linf.csts}}
#'   and optionally refined by \code{\link{refine.linf.csts}} /
#'   \code{\link{refine.linf.csts.iter}}.
#'
#' @return A \code{"linf.csts"} object using the "rare" view.
#'
#' @export
expand.rare <- function(csts) {

  validate.linf.csts(csts)

  if (is.null(csts$cell.label.rare)) {
    stop("expand.rare: no precomputed rare labels found in object")
  }

  if (!is.null(csts$cell.index.rare)) {
    csts$cell.index <- csts$cell.index.rare
  }
  if (!is.null(csts$cell.id.rare)) {
    csts$cell.id <- csts$cell.id.rare
  }
  csts$cell.label <- csts$cell.label.rare
  csts$low.freq.policy <- "rare"

  if (!is.null(csts$cst.id.levels)) {
    if (is.null(csts$cst.id.levels.rare)) {
      depth <- csts$cst.depth
      csts$cst.id.levels[[depth]] <- csts$cell.id %||% csts$cell.label
    } else {
      csts$cst.id.levels <- csts$cst.id.levels.rare
    }
  }

  if (!is.null(csts$cst.levels)) {
    if (is.null(csts$cst.levels.rare)) {
      ## Fall back to updating only the leaf level
      depth <- csts$cst.depth
      csts$cst.levels[[depth]] <- csts$cell.label.rare
    } else {
      csts$cst.levels <- csts$cst.levels.rare
    }
  }

  csts$landmarks <- NULL
  class(csts) <- "linf.csts"
  csts
}

#' Print method for \code{"linf.csts"}
#'
#' @param x A \code{"linf.csts"} object.
#' @param ... Unused.
#' @return The input object, invisibly.
#' @export
print.linf.csts <- function(x, ...) {

  validate.linf.csts(x)

  ## Determine hierarchy
  if (is.null(x$cst.levels)) {
    levels <- list(level1 = x$cell.label)
    max.depth <- 1L
  } else {
    levels <- x$cst.levels
    max.depth <- x$cst.depth
  }

  cat("\n================================================================================\n")
  cat("L-infinity CST Hierarchy\n")
  cat("================================================================================\n")
  cat("Total samples: ", length(x$cell.label), "\n")
  cat("Max depth:     ", max.depth, "\n")

  ## Policy stamp (if available)
  if (!is.null(x$low.freq.policy) && !is.null(x$rare.label)) {
    cat("Low-freq policy: ", x$low.freq.policy,
        " (rare.label: ", x$rare.label, ")\n", sep = "")
  }

  cat("--------------------------------------------------------------------------------\n\n")

  for (d in seq_along(levels)) {
    cat("Depth", d, "\n")
    tab <- sort(table(levels[[d]]), decreasing = TRUE)
    for (nm in names(tab)) {
      cat("  ", nm, ": ", tab[[nm]], "\n", sep = "")
    }
    cat("\n")
  }

  invisible(x)
}

#' Summarize CST hierarchy statistics
#'
#' @param object A \code{"linf.csts"} object.
#' @param ... Unused.
#' @return Data frame with depth-wise statistics
#' @export
summary.linf.csts <- function(object, ...) {

  validate.linf.csts(object)

  ## Determine hierarchy
  if (is.null(object$cst.levels)) {
    levels <- list(level1 = object$cell.label)
  } else {
    levels <- object$cst.levels
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
      median.size = stats::median(tab),
      min.size = min(tab),
      max.size = max(tab)
    ))
  }

  out
}

escape.latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x, fixed = TRUE)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x, perl = TRUE)
  x
}

df.to.latex <- function(df,
                        caption = NULL,
                        label = NULL,
                        digits = 1,
                        include.rownames = FALSE,
                        booktabs = TRUE,
                        use.float = TRUE,
                        print = FALSE) {
  if (!is.data.frame(df)) stop("df.to.latex: df must be a data.frame")

  if (!include.rownames) {
    row_spec <- "l"
  } else {
    row_spec <- "l"
    df <- cbind(Row = rownames(df), df, stringsAsFactors = FALSE)
  }

  col_spec <- paste0(row_spec, paste(rep("r", ncol(df) - 1L), collapse = ""))
  if (ncol(df) == 1L) col_spec <- "l"

  tab_top <- if (booktabs) "\\toprule" else "\\hline"
  tab_mid <- if (booktabs) "\\midrule" else "\\hline"
  tab_bot <- if (booktabs) "\\bottomrule" else "\\hline"

  lines <- character(0)
  if (use.float) lines <- c(lines, "\\begin{table}[htbp]", "\\centering")
  if (!is.null(caption)) lines <- c(lines, paste0("\\caption{", escape.latex(caption), "}"))
  if (!is.null(label)) lines <- c(lines, paste0("\\label{", escape.latex(label), "}"))

  lines <- c(lines, paste0("\\begin{tabular}{", col_spec, "}"), tab_top)

  header <- paste(escape.latex(colnames(df)), collapse = " & ")
  lines <- c(lines, paste0(header, " \\\\"), tab_mid)

  for (i in seq_len(nrow(df))) {
    vals <- vapply(df[i, , drop = FALSE], function(v) {
      if (is.numeric(v)) {
        format(round(v, digits), trim = TRUE, scientific = FALSE)
      } else {
        as.character(v)
      }
    }, character(1))
    vals <- escape.latex(vals)
    lines <- c(lines, paste0(paste(vals, collapse = " & "), " \\\\"))
  }

  lines <- c(lines, tab_bot, "\\end{tabular}")
  if (use.float) lines <- c(lines, "\\end{table}")

  if (isTRUE(print)) cat(paste(lines, collapse = "\n"), "\n")
  lines
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
