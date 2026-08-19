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
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} preserves
#'   sparse input and otherwise uses the dense path.
#'
#' @return Numeric matrix of same dimensions as X, L-infinity normalized.
#'
#' @details
#' Zero rows have undefined L-infinity direction. By convention, they are
#' preserved as all-zero rows and will yield NA labels in downstream
#' dominant-feature or dCST assignment.
#'
#' @export
normalize.linf <- function(X,
                           tol = 0,
                           backend = c("auto", "dense", "sparse")) {
  prep <- linf.prepare.matrix(X, backend = backend, fun.name = "normalize.linf")
  X <- prep$X
  backend <- prep$backend

  linf.validate.matrix(X, backend = backend, fun.name = "normalize.linf")
  if (!is.numeric(tol) || length(tol) != 1L || tol < 0) {
    stop("normalize.linf: tol must be a single non-negative number")
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

linf.normalize.low.freq.policy <- function(low.freq.policy) {
  if (length(low.freq.policy) > 1L) {
    low.freq.policy <- low.freq.policy[[1L]]
  }

  match.arg(low.freq.policy, choices = c("pure", "absorb"))
}

linf.active.low.freq.view <- function(low.freq.policy) {
  if (is.null(low.freq.policy)) {
    return("active")
  }

  if (identical(low.freq.policy, "pure")) {
    return("pure")
  }

  low.freq.policy
}

#' Truncated dominant community state types with configurable low-support handling
#'
#' @description
#' Forms provisional depth-1 dominance sample sets from the dominant feature
#' of each sample and then applies the minimum support threshold \code{n0}.
#' Sets with fewer than \code{n0} samples are handled according to
#' \code{low.freq.policy}:
#' \itemize{
#'   \item \code{"pure"}: retain only sets with support >= \code{n0} as named
#'     dCSTs and collapse all low-support sets into \code{rare.label}.
#'   \item \code{"absorb"}: reassign each low-support sample to the retained
#'     state with the largest value among the retained features (ties handled
#'     by \code{tie.method}).
#' }
#'
#' @param S Numeric matrix (samples x features), typically L-infinity relatives.
#' @param feature.ids Optional character vector of stable feature identifiers,
#'   length \code{ncol(S)}.
#' @param feature.labels Optional character vector of display labels, length
#'   \code{ncol(S)}.
#' @param n0 Integer >= 1. Minimum support required to retain a dominance
#'   sample set.
#' @param low.freq.policy Character. One of \code{"pure"} or \code{"absorb"}.
#'   Default: \code{"pure"}.
#' @param rare.label Character scalar used when \code{low.freq.policy = "pure"}.
#'   Default: \code{"RARE_DOMINANT"}.
#' @param tie.method Character. Tie handling passed to \code{linf.dominant.features()} and used
#'   during absorb reassignment ("first", "random", "error").
#' @param return.diagnostics Logical. If TRUE, return reassignment diagnostics.
#' @param return.landmarks Logical. If TRUE, attach a depth-1 landmark summary
#'   computed by \code{\link{linf.landmarks}}.
#' @param landmark.types Character vector of landmark types passed to
#'   \code{\link{linf.landmarks}} when \code{return.landmarks = TRUE}.
#' @param landmark.view Character. Landmark view passed to
#'   \code{\link{linf.landmarks}} when \code{return.landmarks = TRUE}.
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} preserves
#'   sparse input and otherwise uses the dense path.
#'
#' @return List with:
#'   \itemize{
#'     \item \code{depth1.feature.index}, \code{depth1.feature.id},
#'       \code{depth1.feature.label}: active depth-1 assignment
#'     \item \code{lineage.id}, \code{lineage.label}: active leaf-lineage assignment
#'     \item policy-specific variants of the depth-1 and leaf-lineage fields,
#'       ending in \code{.pure} or \code{.absorb}
#'     \item \code{lineage.ids}, \code{lineage.labels}: active hierarchy, plus
#'       policy-specific \code{.pure} and \code{.absorb} hierarchies
#'     \item \code{depth}: current hierarchy depth
#'     \item \code{retained.feature.indices}, \code{retained.feature.ids}, \code{retained.feature.labels}
#'     \item \code{provisional.feature.index}, \code{provisional.feature.id}, \code{provisional.feature.label}
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
                      low.freq.policy = c("pure", "absorb"),
                      rare.label = "RARE_DOMINANT",
                      tie.method = c("first", "random", "error"),
                      return.diagnostics = FALSE,
                      return.landmarks = FALSE,
                      landmark.types = c("endpoint.max", "endpoint.min"),
                      landmark.view = c("active", "pure", "absorb"),
                      backend = c("auto", "dense", "sparse")) {

    low.freq.policy <- linf.normalize.low.freq.policy(low.freq.policy)
    tie.method <- match.arg(tie.method)
    landmark.view <- match.arg(landmark.view)
    prep <- linf.prepare.matrix(S, backend = backend, fun.name = "linf.csts")
    X <- prep$X
    backend <- prep$backend

    if (!is.numeric(n0) || length(n0) != 1L || n0 < 1 || n0 %% 1 != 0) {
        stop("linf.csts: n0 must be integer >= 1")
    }
    if (!is.character(rare.label) || length(rare.label) != 1L || !nzchar(rare.label)) {
        stop("linf.csts: rare.label must be a non-empty character scalar")
    }

    linf.validate.matrix(X, backend = backend, fun.name = "linf.csts")

    meta <- resolve.linf.feature.meta(X, feature.ids = feature.ids, feature.labels = feature.labels)
    fid <- meta$feature.ids
    lev <- meta$feature.labels

    raw <- linf.dominant.features(X,
                      feature.ids = fid,
                      feature.labels = lev,
                      tie.method = tie.method,
                      backend = backend)
    tab <- sort(table(raw$label[!is.na(raw$label)]), decreasing = TRUE)
    tab.id <- sort(table(raw$id[!is.na(raw$id)]), decreasing = TRUE)

    kept.lbl <- names(tab[tab >= n0])
    kept.idx <- match(kept.lbl, lev)
    kept.id <- fid[kept.idx]

    n <- nrow(X)

    ## Pure-policy labels: retain only states with support >= n0.
    is.kept <- !is.na(raw$label) & (raw$label %in% kept.lbl)

    depth1.feature.idx.pure <- raw$index
    lineage.id.pure <- raw$id
    lineage.label.pure <- raw$label
    depth1.feature.idx.pure[!is.kept] <- NA_integer_
    lineage.id.pure[!is.kept] <- rare.label
    lineage.label.pure[!is.kept] <- rare.label
    depth1.feature.id.pure <- ifelse(
        is.na(depth1.feature.idx.pure), NA_character_, fid[depth1.feature.idx.pure]
    )
    depth1.feature.label.pure <- ifelse(
        is.na(depth1.feature.idx.pure), NA_character_, lev[depth1.feature.idx.pure]
    )

    ## Absorb-policy labels: reassign low-support (and zero-row) samples
    ## into retained states.
    depth1.feature.idx.absorb <- raw$index
    lineage.id.absorb <- raw$id
    lineage.label.absorb <- raw$label

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
            if (backend == "sparse") {
                new.idx <- linf.absorb.sparse(
                    X,
                    to.absorb = to.absorb,
                    kept.idx = kept.idx,
                    tie.method = tie.method,
                    fallback.idx = fallback.idx
                )
            } else {
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
            }

            reassigned[to.absorb] <- TRUE
            reassigned.from[to.absorb] <- raw$label[to.absorb]
            reassigned.to[to.absorb]   <- lev[new.idx]

            depth1.feature.idx.absorb[to.absorb] <- new.idx
            lineage.id.absorb[to.absorb] <- fid[new.idx]
            lineage.label.absorb[to.absorb] <- lev[new.idx]
        }

    } else {

        ## No retained states at this n0:
        ## - pure policy: everyone is rare.label (already set above)
        ## - absorb-policy: undefined; keep NA labels
        depth1.feature.idx.absorb[] <- NA_integer_
        lineage.id.absorb[] <- NA_character_
        lineage.label.absorb[] <- NA_character_
    }

    ## Select active labeling
    if (low.freq.policy == "pure") {
        depth1.feature.idx <- depth1.feature.idx.pure
        depth1.feature.id <- depth1.feature.id.pure
        depth1.feature.label <- depth1.feature.label.pure
        lineage.id <- lineage.id.pure
        lineage.label <- lineage.label.pure
    } else {
        depth1.feature.idx <- depth1.feature.idx.absorb
        depth1.feature.id <- lineage.id.absorb
        depth1.feature.label <- lineage.label.absorb
        lineage.id <- lineage.id.absorb
        lineage.label <- lineage.label.absorb
    }

    names(depth1.feature.idx) <- rownames(X)
    names(depth1.feature.id) <- rownames(X)
    names(depth1.feature.label) <- rownames(X)
    names(lineage.id) <- rownames(X)
    names(lineage.label) <- rownames(X)

    names(depth1.feature.idx.pure) <- rownames(X)
    names(depth1.feature.id.pure) <- rownames(X)
    names(depth1.feature.label.pure) <- rownames(X)
    names(lineage.id.pure) <- rownames(X)
    names(lineage.label.pure) <- rownames(X)

    names(depth1.feature.idx.absorb) <- rownames(X)
    names(lineage.id.absorb) <- rownames(X)
    names(lineage.label.absorb) <- rownames(X)

    out <- list(
        depth1.feature.index       = depth1.feature.idx,
        depth1.feature.id          = depth1.feature.id,
        depth1.feature.label       = depth1.feature.label,
        lineage.id          = lineage.id,
        lineage.label       = lineage.label,
        depth1.feature.index.pure  = depth1.feature.idx.pure,
        depth1.feature.id.pure     = depth1.feature.id.pure,
        depth1.feature.label.pure  = depth1.feature.label.pure,
        lineage.id.pure     = lineage.id.pure,
        lineage.label.pure  = lineage.label.pure,
        depth1.feature.index.absorb = depth1.feature.idx.absorb,
        depth1.feature.id.absorb    = lineage.id.absorb,
        depth1.feature.label.absorb = lineage.label.absorb,
        lineage.id.absorb   = lineage.id.absorb,
        lineage.label.absorb = lineage.label.absorb,
        retained.feature.indices   = kept.idx,
        retained.feature.ids    = kept.id,
        retained.feature.labels   = kept.lbl,
        provisional.feature.index        = raw$index,
        provisional.feature.id           = raw$id,
        provisional.feature.label        = raw$label,
        size.table       = tab,
        size.table.id    = tab.id,
        feature.ids      = fid,
        feature.labels   = lev,
        matrix.backend   = backend,
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

    out$lineage.ids <- list(level1 = out$lineage.id)
    out$lineage.ids.pure <- list(level1 = out$lineage.id.pure)
    out$lineage.ids.absorb <- list(level1 = out$lineage.id.absorb)
    out$lineage.labels <- list(level1 = out$lineage.label)
    out$lineage.labels.pure <- list(level1 = out$lineage.label.pure)
    out$lineage.labels.absorb <- list(level1 = out$lineage.label.absorb)
    out$depth <- 1L

    class(out) <- "linf.csts"

    if (isTRUE(return.landmarks)) {
        out$landmarks <- linf.landmarks(
            X,
            out,
            depth = 1L,
            view = landmark.view,
            landmark.types = landmark.types,
            tie.method = tie.method,
            backend = backend
        )
    }

    out
}

#' From counts to filtered, normalized, truncated dCSTs
#'
#' Convenience pipeline: runs \code{filter.asv()} on counts, computes
#' L-infinity-normalized relatives on the filtered counts, and assigns
#' truncated dCSTs.
#'
#' @param S.counts Numeric matrix of counts (samples x features).
#' @param ... Arguments passed to \code{filter.asv()} (e.g., \code{min.lib},
#'   \code{prev.prop}, \code{min.count}, \code{min.rel}).
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}.
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
asv.to.linf.csts <- function(S.counts,
                             ...,
                             backend = c("auto", "dense", "sparse")) {
  filt <- filter.asv(S.counts, ...)
  M <- normalize.linf(filt$counts, backend = backend)
  cst <- linf.csts(M, backend = backend)
  list(filter = filt, linf.rel = M, csts = cst)
}

validate.linf.csts <- function(obj) {

  stopifnot(is.list(obj))
  stopifnot(!is.null(obj$lineage.label))

  n <- length(obj$lineage.label)

  ## Case 1: explicit hierarchy (refined dCSTs)
  if (!is.null(obj$lineage.labels)) {

    stopifnot(is.list(obj$lineage.labels))
    stopifnot(!is.null(obj$depth))
    stopifnot(obj$depth == length(obj$lineage.labels))

    for (lvl in obj$lineage.labels) {
      stopifnot(length(lvl) == n)
    }

    return(invisible(TRUE))
  }

  ## A depth-1 dCST has no explicit hierarchy.
  return(invisible(TRUE))
}

#' Refine a dCST hierarchy by one level
#'
#' @description
#' Selects well-supported leaf dominance-lineages and refines them by dropping
#' the dominant feature(s) encoded in the parent label path and re-applying
#' \code{\link{linf.csts}} to the remaining features. The resulting child
#' labels are appended to the parent label using \code{sep}.
#'
#' Low-support child lineages are handled by \code{low.freq.policy}. When
#' \code{low.freq.policy = "pure"}, rare buckets at depth >= 2 become
#' parent-prefixed automatically via the hierarchical \code{paste(parent, child, sep = sep)}.
#'
#' @param M Numeric matrix (samples x features) used for refinement. Column names
#'   should match feature labels used in dCST names.
#' @param csts A \code{"linf.csts"} object.
#' @param n0 Integer >= 1. Minimum support required to retain a child lineage
#'   (passed to \code{linf.csts}).
#' @param refinement.factor Numeric > 0. Auto-refine parent lineages with
#'   support >= \code{refinement.factor * n0}.
#' @param sep Character scalar used to concatenate hierarchical labels.
#' @param low.freq.policy Character. One of \code{"pure"} or \code{"absorb"}.
#'   Default: \code{"pure"}.
#' @param rare.label Character scalar for rare buckets when \code{low.freq.policy = "pure"}.
#' @param verbose Logical. If TRUE, print progress information.
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} inherits the
#'   backend from \code{M} or from \code{csts} when available.
#'
#' @return Updated \code{"linf.csts"} object with \code{depth} increased by one and
#'   updated \code{lineage.label}. Policy-specific views are stored in
#'   \code{lineage.label.pure} and \code{lineage.label.absorb}.
#'
#' @export
refine.linf.csts <- function(M,
                             csts,
                             n0 = 50,
                             refinement.factor = 2,
                             sep = "__",
                             low.freq.policy = c("pure", "absorb"),
                             rare.label = "RARE_DOMINANT",
                             verbose = TRUE,
                             backend = c("auto", "dense", "sparse")) {

    validate.linf.csts(csts)

    low.freq.policy <- linf.normalize.low.freq.policy(low.freq.policy)
    if (missing(backend) && !is.null(csts$matrix.backend)) {
        backend <- csts$matrix.backend
    }
    prep <- linf.prepare.matrix(M, backend = backend, fun.name = "refine.linf.csts")
    M <- prep$X
    backend <- prep$backend
    linf.validate.matrix(M, backend = backend, fun.name = "refine.linf.csts")

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
    if (is.null(csts$lineage.labels)) {
        csts$lineage.labels <- list(level1 = csts$lineage.label)
        csts$depth <- 1L
        csts$sep <- sep

        ## Ensure both policy-specific views exist at depth 1
        if (is.null(csts$lineage.label.pure)) {
            csts$lineage.label.pure <- csts$lineage.label
        }
        if (is.null(csts$lineage.label.absorb)) {
            csts$lineage.label.absorb <- csts$lineage.label
        }

        csts$lineage.labels.pure <- list(level1 = csts$lineage.label.pure)
        csts$lineage.labels.absorb <- list(level1 = csts$lineage.label.absorb)
    }

    if (is.null(csts$lineage.ids)) {
        csts$lineage.ids <- list(level1 = csts$lineage.id %||% csts$lineage.label)
    }

    if (is.null(csts$feature.ids)) {
        csts$feature.ids <- colnames(M)
    }
    if (is.null(csts$feature.labels)) {
        csts$feature.labels <- csts$feature.ids
    }

    ## Ensure both stored policy views are present.
    if (is.null(csts$lineage.labels.pure)) csts$lineage.labels.pure <- csts$lineage.labels
    if (is.null(csts$lineage.labels.absorb)) csts$lineage.labels.absorb <- csts$lineage.labels
    if (is.null(csts$lineage.ids.pure)) csts$lineage.ids.pure <- csts$lineage.ids
    if (is.null(csts$lineage.ids.absorb)) csts$lineage.ids.absorb <- csts$lineage.ids

    depth <- csts$depth + 1L
    parent.ids <- csts$lineage.ids[[depth - 1L]]
    parent.ids.pure <- csts$lineage.ids.pure[[depth - 1L]]
    parent.ids.absorb <- csts$lineage.ids.absorb[[depth - 1L]]
    parent.labels <- csts$lineage.labels[[depth - 1L]]

    parent.labels.pure <- csts$lineage.labels.pure[[depth - 1L]]
    parent.labels.absorb <- csts$lineage.labels.absorb[[depth - 1L]]

    refined.ids.pure <- parent.ids.pure
    refined.ids.absorb <- parent.ids.absorb
    refined.labels.pure <- parent.labels.pure
    refined.labels.absorb <- parent.labels.absorb

    lineage.sizes <- sort(table(parent.ids), decreasing = TRUE)

    threshold <- refinement.factor * n0
    refine.lineages <- names(lineage.sizes[lineage.sizes >= threshold])

    if (verbose) {
        cat("========================================\n")
        cat("AUTO-REFINEMENT MODE\n")
        cat("========================================\n")
        cat("Refinement threshold:", threshold, "\n")
        cat("Dominance-lineages selected for refinement:", length(refine.lineages), "\n\n")
    }

    for (lineage in refine.lineages) {
        ## Rare buckets are synthetic catch-all groups, not real feature paths.
        ## They should never be refined further.
        if (identical(lineage, rare.label)) next

        idx <- which(parent.ids == lineage)
        parent.id.pure <- parent.ids.pure[idx[1L]]
        parent.id.absorb <- parent.ids.absorb[idx[1L]]
        parent.label.pure <- parent.labels.pure[idx[1L]]
        parent.label.absorb <- parent.labels.absorb[idx[1L]]

        parent.features <- strsplit(lineage, sep, fixed = TRUE)[[1]]
        drop.idx <- match(parent.features, csts$feature.ids)
        drop.idx <- drop.idx[!is.na(drop.idx)]

        ## If no parent features match columns (e.g., the lineage is a rare
        ## bucket), do not drop any columns.
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
                             rare.label = rare.label,
                             backend = backend)

        sub.ids.pure <- sub.csts$lineage.id.pure %||% sub.csts$lineage.label.pure
        sub.ids.absorb <- sub.csts$lineage.id.absorb %||% sub.csts$lineage.label.absorb
        sub.labels.pure <- sub.csts$lineage.label.pure
        sub.labels.absorb <- sub.csts$lineage.label.absorb

        ## If refinement yields no retained child lineages, skip.
        if (all(sub.labels.pure == rare.label)) next

        refined.ids.pure[idx] <- paste(parent.id.pure, sub.ids.pure, sep = sep)
        refined.ids.absorb[idx] <- paste(parent.id.absorb, sub.ids.absorb, sep = sep)
        refined.labels.pure[idx] <- paste(parent.label.pure, sub.labels.pure, sep = sep)
        refined.labels.absorb[idx] <- paste(parent.label.absorb, sub.labels.absorb, sep = sep)
    }

    ## Update dCST object (store both views; active view chosen by low.freq.policy)
    csts$lineage.ids[[depth]] <- if (low.freq.policy == "pure") refined.ids.pure else refined.ids.absorb
    csts$lineage.ids.pure[[depth]] <- refined.ids.pure
    csts$lineage.ids.absorb[[depth]] <- refined.ids.absorb
    csts$lineage.labels.pure[[depth]] <- refined.labels.pure
    csts$lineage.labels.absorb[[depth]] <- refined.labels.absorb

    if (low.freq.policy == "pure") {
        refined.labels <- refined.labels.pure
    } else {
        refined.labels <- refined.labels.absorb
    }

    csts$lineage.labels[[depth]] <- refined.labels
    csts$depth <- depth
    csts$sep <- sep

    csts$lineage.id.pure <- refined.ids.pure
    csts$lineage.id.absorb <- refined.ids.absorb
    csts$lineage.id <- csts$lineage.ids[[depth]]
    csts$lineage.label.pure <- refined.labels.pure
    csts$lineage.label.absorb <- refined.labels.absorb
    csts$lineage.label <- refined.labels

    csts$low.freq.policy <- low.freq.policy
    csts$rare.label <- rare.label
    csts$matrix.backend <- backend
    csts$landmarks <- NULL

    class(csts) <- "linf.csts"

    csts
}

#' Iteratively refine dCSTs by one additional level
#'
#' @description
#' Appends one refinement level to an existing dCST hierarchy produced by
#' \code{\link{refine.linf.csts}} or \code{\link{refine.linf.csts.iter}}.
#' Dominance-lineages to refine can be provided explicitly or selected
#' automatically based on \code{refinement.factor * n0}.
#'
#' @param M Numeric matrix used for refinement.
#' @param refined A \code{"linf.csts"} object with an existing hierarchy.
#' @param lineages.to.refine Character vector of leaf dominance-lineage IDs to
#'   refine, or NULL for auto-selection.
#' @param n0 Integer >= 1. Minimum support required to retain a child lineage.
#' @param refinement.factor Numeric > 0. Auto-selection threshold multiplier.
#' @param sep Character scalar used to concatenate hierarchical labels.
#' @param low.freq.policy Character. One of \code{"pure"} or \code{"absorb"}.
#'   Default: \code{"pure"}.
#' @param rare.label Character scalar for rare buckets when \code{low.freq.policy = "pure"}.
#' @param verbose Logical. If TRUE, print progress information.
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} inherits the
#'   backend from \code{M} or from \code{refined} when available.
#'
#' @return Updated \code{"linf.csts"} object with \code{depth} increased by one.
#'
#' @export
refine.linf.csts.iter <- function(M,
                                  refined,
                                  lineages.to.refine = NULL,
                                  n0 = 50,
                                  refinement.factor = 5,
                                  sep = "__",
                                  low.freq.policy = c("pure", "absorb"),
                                  rare.label = "RARE_DOMINANT",
                                  verbose = TRUE,
                                  backend = c("auto", "dense", "sparse")) {

    validate.linf.csts(refined)

    low.freq.policy <- linf.normalize.low.freq.policy(low.freq.policy)
    if (missing(backend) && !is.null(refined$matrix.backend)) {
        backend <- refined$matrix.backend
    }
    prep <- linf.prepare.matrix(M, backend = backend, fun.name = "refine.linf.csts.iter")
    M <- prep$X
    backend <- prep$backend
    linf.validate.matrix(M, backend = backend, fun.name = "refine.linf.csts.iter")

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

    ## Ensure both stored policy views are present.
    if (is.null(refined$lineage.labels.pure)) refined$lineage.labels.pure <- refined$lineage.labels
    if (is.null(refined$lineage.labels.absorb)) refined$lineage.labels.absorb <- refined$lineage.labels
    if (is.null(refined$lineage.ids)) refined$lineage.ids <- list(level1 = refined$lineage.id %||% refined$lineage.label)
    if (is.null(refined$lineage.ids.pure)) refined$lineage.ids.pure <- refined$lineage.ids
    if (is.null(refined$lineage.ids.absorb)) refined$lineage.ids.absorb <- refined$lineage.ids
    if (is.null(refined$feature.ids)) refined$feature.ids <- colnames(M)
    if (is.null(refined$feature.labels)) refined$feature.labels <- refined$feature.ids

    depth <- refined$depth

    parent.ids <- refined$lineage.ids[[depth]]
    parent.ids.pure <- refined$lineage.ids.pure[[depth]]
    parent.ids.absorb <- refined$lineage.ids.absorb[[depth]]
    parent.labels <- refined$lineage.labels[[depth]]
    parent.labels.pure <- refined$lineage.labels.pure[[depth]]
    parent.labels.absorb <- refined$lineage.labels.absorb[[depth]]

    lineage.sizes <- sort(table(parent.ids), decreasing = TRUE)
    threshold <- refinement.factor * n0

    if (is.null(lineages.to.refine)) {
        lineages.to.refine <- names(lineage.sizes[lineage.sizes >= threshold])
    }

    if (verbose) {
        cat("Auto-selecting dominance-lineages for iterative refinement (threshold =", threshold, "):\n")
        for (c in lineages.to.refine) {
            cat(" -", c, ":", lineage.sizes[c], "samples\n")
        }
        cat("\n")
    }

    new.ids.pure <- parent.ids.pure
    new.ids.absorb <- parent.ids.absorb
    new.labels.pure <- parent.labels.pure
    new.labels.absorb <- parent.labels.absorb

    for (lineage in lineages.to.refine) {
        ## Rare buckets are synthetic catch-all groups, not real feature paths.
        ## They should never be refined further.
        if (identical(lineage, rare.label)) next

        idx <- which(parent.ids == lineage)
        parent.id.pure <- parent.ids.pure[idx[1L]]
        parent.id.absorb <- parent.ids.absorb[idx[1L]]
        parent.label.pure <- parent.labels.pure[idx[1L]]
        parent.label.absorb <- parent.labels.absorb[idx[1L]]

        parents <- strsplit(lineage, sep, fixed = TRUE)[[1]]
        drop.idx <- match(parents, refined$feature.ids)
        drop.idx <- drop.idx[!is.na(drop.idx)]

        ## If no parent features match columns (e.g., the lineage is a rare
        ## bucket), do not drop any columns.
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
                             rare.label = rare.label,
                             backend = backend)

        sub.ids.pure <- sub.csts$lineage.id.pure %||% sub.csts$lineage.label.pure
        sub.ids.absorb <- sub.csts$lineage.id.absorb %||% sub.csts$lineage.label.absorb
        sub.labels.pure <- sub.csts$lineage.label.pure
        sub.labels.absorb <- sub.csts$lineage.label.absorb

        ## If refinement yields no retained child lineages, skip.
        if (all(sub.labels.pure == rare.label)) next

        new.ids.pure[idx] <- paste(parent.id.pure, sub.ids.pure, sep = sep)
        new.ids.absorb[idx] <- paste(parent.id.absorb, sub.ids.absorb, sep = sep)
        new.labels.pure[idx] <- paste(parent.label.pure, sub.labels.pure, sep = sep)
        new.labels.absorb[idx] <- paste(parent.label.absorb, sub.labels.absorb, sep = sep)
    }

    refined$lineage.ids[[depth + 1L]] <- if (low.freq.policy == "pure") new.ids.pure else new.ids.absorb
    refined$lineage.ids.pure[[depth + 1L]] <- new.ids.pure
    refined$lineage.ids.absorb[[depth + 1L]] <- new.ids.absorb
    refined$lineage.labels.pure[[depth + 1L]] <- new.labels.pure
    refined$lineage.labels.absorb[[depth + 1L]] <- new.labels.absorb

    if (low.freq.policy == "pure") {
        new.labels <- new.labels.pure
    } else {
        new.labels <- new.labels.absorb
    }

    refined$lineage.labels[[depth + 1L]] <- new.labels
    refined$depth <- depth + 1L

    refined$lineage.id.pure <- new.ids.pure
    refined$lineage.id.absorb <- new.ids.absorb
    refined$lineage.id <- refined$lineage.ids[[depth + 1L]]
    refined$lineage.label.pure <- new.labels.pure
    refined$lineage.label.absorb <- new.labels.absorb
    refined$lineage.label <- new.labels

    refined$low.freq.policy <- low.freq.policy
    refined$rare.label <- rare.label
    refined$matrix.backend <- backend
    refined$sep <- sep
    refined$landmarks <- NULL

    class(refined) <- "linf.csts"

    refined
}

#' Select a stored dCST policy view
#'
#' @description
#' Returns a copy of a \code{"linf.csts"} object using either the stored
#' \code{"pure"} or \code{"absorb"} hierarchy. This does not recompute dCSTs.
#'
#' @param csts A \code{"linf.csts"} object produced by \code{\link{linf.csts}}
#'   and optionally refined by \code{\link{refine.linf.csts}} /
#'   \code{\link{refine.linf.csts.iter}}.
#' @param view Character. The stored policy view to activate: \code{"pure"} or
#'   \code{"absorb"}.
#'
#' @return A \code{"linf.csts"} object using the requested view.
#'
#' @export
dcst.view <- function(csts, view = c("absorb", "pure")) {

  validate.linf.csts(csts)
  view <- match.arg(view)

  lineage.label <- csts[[paste0("lineage.label.", view)]]
  lineage.id <- csts[[paste0("lineage.id.", view)]]
  depth1.feature.index <- csts[[paste0("depth1.feature.index.", view)]]
  depth1.feature.id <- csts[[paste0("depth1.feature.id.", view)]]
  depth1.feature.label <- csts[[paste0("depth1.feature.label.", view)]]
  lineage.ids <- csts[[paste0("lineage.ids.", view)]]
  lineage.labels <- csts[[paste0("lineage.labels.", view)]]

  if (is.null(lineage.label) || is.null(lineage.ids) || is.null(lineage.labels)) {
    stop(sprintf('dcst.view: object does not contain the "%s" view', view))
  }

  csts$depth1.feature.index <- depth1.feature.index
  csts$depth1.feature.id <- depth1.feature.id
  csts$depth1.feature.label <- depth1.feature.label
  csts$lineage.id <- lineage.id
  csts$lineage.label <- lineage.label
  csts$lineage.ids <- lineage.ids
  csts$lineage.labels <- lineage.labels
  csts$low.freq.policy <- view

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
  if (is.null(x$lineage.labels)) {
    levels <- list(level1 = x$lineage.label)
    max.depth <- 1L
  } else {
    levels <- x$lineage.labels
    max.depth <- x$depth
  }

  cat("\n================================================================================\n")
  cat("Dominant Community State Type Hierarchy\n")
  cat("================================================================================\n")
  cat("Total samples: ", length(x$lineage.label), "\n")
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

#' Summarize dCST hierarchy statistics
#'
#' @param object A \code{"linf.csts"} object.
#' @param ... Unused.
#' @return Data frame with depth-wise statistics
#' @export
summary.linf.csts <- function(object, ...) {

  validate.linf.csts(object)

  ## Determine hierarchy
  if (is.null(object$lineage.labels)) {
    levels <- list(level1 = object$lineage.label)
  } else {
    levels <- object$lineage.labels
  }

  out <- data.frame(
    depth = integer(0),
    n.lineages = integer(0),
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
      n.lineages = length(tab),
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

#' Render dCSTs as a LaTeX table
#'
#' @description
#' Generates LaTeX code summarizing dCSTs at a specified hierarchy depth.
#' Uses the explicit dCST hierarchy stored in \code{lineage.labels} when available.
#'
#' @param csts A dCST object produced by \code{linf.csts()} and optionally refined
#'   by \code{refine.linf.csts()} or \code{refine.linf.csts.iter()}.
#' @param depth Integer. dCST depth to render. Default is the leaf level
#'   (\code{csts$depth}). If \code{lineage.labels} is missing, depth is ignored.
#' @param caption Character or NULL. LaTeX table caption.
#' @param label Character or NULL. LaTeX label for referencing the table.
#' @param digits Integer. Digits to display for percentages. Default: 1.
#' @param include.percent Logical. If TRUE, include percent of total samples.
#'
#' @return Character vector containing LaTeX table code.
#'
#' @details
#' This function is hierarchy-aware. It does not infer dCST depth from label
#' strings. If the dCST object does not contain \code{lineage.labels}, the function
#' falls back to using \code{lineage.label} as a single-level dCST.
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
  if (!is.null(csts$lineage.labels)) {
    if (is.null(depth)) {
      depth <- csts$depth
    }
    if (depth < 1 || depth > csts$depth) {
      stop("latex.linf.csts: invalid depth")
    }
    labels <- csts$lineage.labels[[depth]]
  } else {
    ## Backward compatibility: flat dCST
    labels <- csts$lineage.label
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
