#' Refine a dCST hierarchy by one level
#'
#' @description
#' Selects leaf dominance-lineages and refines them by dropping the dominant
#' feature(s) stored in explicit node paths and re-applying
#' \code{\link{linf.csts}} to the remaining features. The resulting child
#' labels are joined using \code{sep}; colliding readable paths receive a node suffix.
#'
#' Low-support child lineages are handled by \code{low.freq.policy}. When
#' \code{low.freq.policy = "pure"}, rare buckets at depth >= 2 become
#' parent-prefixed automatically via the hierarchical \code{paste(parent, child, sep = sep)}.
#'
#' @param M Numeric matrix (samples x features) used for refinement. Columns
#'   must correspond, in order, to the stable feature IDs stored in \code{csts}.
#'   Rows must remain in fitted sample order. Both dimensions are checked;
#'   named axes are checked against `csts$input.dimnames` when available.
#'   Unnamed axes/legacy fits use positional matching. No automatic reordering.
#' @param csts A \code{"linf.csts"} object.
#' @param lineages.to.refine Optional character vector of leaf
#'   dominance-lineage IDs to refine. When \code{NULL}, lineages are selected
#'   automatically using \code{refinement.factor * n0}.
#' @param n0 Integer >= 1. Minimum support required to retain a child lineage
#'   (passed to \code{linf.csts}).
#' @param refinement.factor Numeric > 0. Auto-refine parent lineages with
#'   support >= \code{refinement.factor * n0}.
#' @param sep Character scalar used to concatenate hierarchical labels.
#' @param low.freq.policy Character. One of \code{"pure"} or \code{"absorb"}.
#'   Default: \code{"pure"}.
#' @param rare.label Character scalar for rare buckets when \code{low.freq.policy = "pure"}.
#' @param verbose Logical. If TRUE, emit progress messages.
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} inherits the
#'   backend from \code{M} or from \code{csts} when available.
#'
#' @return Updated \code{"linf.csts"} object with \code{depth} increased by one and
#'   updated \code{lineage.label} and a new `history` entry recording the
#'   settings and selected/refined nodes at this depth. Policy-specific views are stored in
#'   \code{lineage.label.pure} and \code{lineage.label.absorb}.
#'
#' @details
#' Selection and child fitting use active leaf membership. Explicit selection
#' bypasses the automatic parent-support rule but not the child threshold.
#' Every call appends a stored level; unselected parents and parents without
#' retained children keep their labels there. Child fits use first-column tie
#' handling; the parent fit's random-tie setting is not inherited. Keep policy,
#' rare label and separator consistent across repeated refinements. Stored
#' synthetic or undefined parents remain terminal in their own policy view,
#' even when the other view is refined. A real feature whose ID equals the
#' rare label is not a synthetic parent and can be refined.
#'
#' @examples
#' M <- rbind(
#'   s1 = c(A = 1.0, B = 0.8, C = 0.1),
#'   s2 = c(A = 1.0, B = 0.7, C = 0.2),
#'   s3 = c(A = 1.0, B = 0.6, C = 0.3),
#'   s4 = c(A = 0.2, B = 1.0, C = 0.8),
#'   s5 = c(A = 0.1, B = 1.0, C = 0.7),
#'   s6 = c(A = 0.3, B = 1.0, C = 0.6)
#' )
#' depth1 <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")
#' depth2 <- refine.linf.csts(
#'   M,
#'   depth1,
#'   lineages.to.refine = "A",
#'   n0 = 2,
#'   low.freq.policy = "absorb",
#'   verbose = FALSE
#' )
#' depth2$lineage.labels[[2]]
#'
#' @export
refine.linf.csts <- function(M,
                             csts,
                             lineages.to.refine = NULL,
                             n0 = 50,
                             refinement.factor = 2,
                             sep = "__",
                             low.freq.policy = c("pure", "absorb"),
                             rare.label = "RARE_DOMINANT",
                             verbose = TRUE,
                             backend = c("auto", "dense", "sparse")) {

    validate.linf.csts(csts)
    csts <- linf.ensure.nodes(csts)

    low.freq.policy <- linf.normalize.low.freq.policy(low.freq.policy)
    if (missing(backend)) {
        backend <- csts$matrix.backend
    }
    prep <- linf.prepare.matrix(M, backend = backend, fun.name = "refine.linf.csts")
    M <- prep$X
    backend <- prep$backend
    linf.validate.matrix(M, backend = backend, fun.name = "refine.linf.csts")
    linf.validate.fit.matrix(M, csts, "refine.linf.csts")

    if (!is.numeric(n0) || length(n0) != 1L || !is.finite(n0) || n0 < 1 || n0 %% 1 != 0) {
        stop("refine.linf.csts: n0 must be integer >= 1")
    }
    if (!is.numeric(refinement.factor) || length(refinement.factor) != 1L ||
        !is.finite(refinement.factor) || refinement.factor <= 0) {
        stop("refine.linf.csts: refinement.factor must be a finite numeric > 0")
    }
    if (!is.character(sep) || length(sep) != 1L || is.na(sep) || !nzchar(sep)) {
        stop("refine.linf.csts: sep must be a non-empty character scalar")
    }
    if (!is.character(rare.label) || length(rare.label) != 1L || is.na(rare.label) || !nzchar(rare.label)) {
        stop("refine.linf.csts: rare.label must be a non-empty character scalar")
    }

    depth <- csts$depth + 1L
    parent.ids <- csts$lineage.id
    active.nodes <- csts$nodes[[depth - 1L]]
    active.paths <- linf.node.paths(csts, csts$low.freq.policy, depth - 1L)
    paths.pure <- linf.node.paths(csts, "pure", depth - 1L)
    paths.absorb <- linf.node.paths(csts, "absorb", depth - 1L)
    rare.ids <- active.nodes$lineage.id[active.nodes$is.rare]

    lineage.sizes <- sort(table(parent.ids), decreasing = TRUE)

    threshold <- refinement.factor * n0
    if (is.null(lineages.to.refine)) {
        lineages.to.refine <- names(lineage.sizes[lineage.sizes >= threshold])
        lineages.to.refine <- setdiff(lineages.to.refine, rare.ids)
        selection.mode <- "automatic"
    } else {
        if (!is.character(lineages.to.refine) || anyNA(lineages.to.refine) ||
            any(!nzchar(lineages.to.refine))) {
            stop("refine.linf.csts: lineages.to.refine must be NULL or a character vector of lineage IDs")
        }
        lineages.to.refine <- unique(lineages.to.refine)
        unknown <- setdiff(lineages.to.refine, names(lineage.sizes))
        if (length(unknown)) {
            stop(
                "refine.linf.csts: unknown lineage ID(s): ",
                paste(unknown, collapse = ", ")
            )
        }
        if (any(lineages.to.refine %in% rare.ids)) {
            stop("refine.linf.csts: the synthetic rare lineage cannot be refined")
        }
        selection.mode <- "explicit"
    }

    if (verbose) {
        progress <- c(
            "========================================",
            paste(toupper(selection.mode), "REFINEMENT MODE"),
            "========================================"
        )
        if (selection.mode == "automatic") {
            progress <- c(progress, paste("Refinement threshold:", threshold))
        }
        progress <- c(
            progress,
            paste(
                "Dominance-lineages selected for refinement:",
                length(lineages.to.refine)
            )
        )
        message(paste(progress, collapse = "\n"))
    }

    refined.node.ids <- character()
    for (lineage in lineages.to.refine) {
        idx <- which(parent.ids == lineage)
        drop.idx <- active.paths[[idx[1L]]]
        drop.idx <- drop.idx[drop.idx > 0L]
        keep.idx <- setdiff(seq_len(ncol(M)), drop.idx)

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
                             feature.ids = csts$feature.ids[keep.idx],
                             feature.labels = csts$feature.labels[keep.idx],
                             n0 = n0,
                             low.freq.policy = low.freq.policy,
                             rare.label = rare.label,
                             backend = backend)

        if (!length(sub.csts$retained.feature.indices)) next
        refined.node.ids <- c(refined.node.ids, csts$node.ids[[depth - 1L]][idx[1L]])
        for (j in seq_along(idx)) {
            row <- idx[[j]]
            # A stored synthetic/undefined parent remains terminal in its view.
            for (view in c("pure", "absorb")) {
                paths <- if (view == "pure") paths.pure else paths.absorb
                parent <- paths[[row]]
                if (!length(parent) || utils::tail(parent, 1L) == 0L) next
                child <- sub.csts[[paste0("depth1.feature.index.", view)]][[j]]
                token <- if (is.na(child)) 0L else as.integer(keep.idx[[child]])
                paths[[row]] <- c(parent, token)
                if (view == "pure") paths.pure <- paths else paths.absorb <- paths
            }
        }
    }

    csts$history <- linf.fit.history(csts)
    csts$history[[depth]] <- list(
        depth = depth, available = TRUE, operation = "refine", n0 = as.integer(n0),
        refinement.factor = refinement.factor, selection.mode = selection.mode,
        selected.lineage.ids = lineages.to.refine,
        selected.node.ids = active.nodes$node.id[match(lineages.to.refine, active.nodes$lineage.id)],
        refined.node.ids = unname(refined.node.ids),
        source.view = csts$low.freq.policy, low.freq.policy = low.freq.policy,
        rare.label = rare.label, sep = sep, tie.method = "first", backend = backend
    )
    csts$depth <- depth
    csts$sep <- sep
    csts$low.freq.policy <- low.freq.policy
    csts$rare.label <- rare.label
    csts <- linf.store.node.level(csts, "pure", paths.pure, depth, sep)
    csts <- linf.store.node.level(csts, "absorb", paths.absorb, depth, sep)
    csts <- linf.activate.nodes(csts)
    csts$matrix.backend <- backend
    csts$landmarks <- NULL

    class(csts) <- "linf.csts"

    csts
}
