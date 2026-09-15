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
#'   Positive ties use original matrix column order for \code{"first"} on both
#'   backends. If all retained values are zero, absorption uses the retained
#'   state with greatest reference support (then lexical label order), regardless
#'   of \code{tie.method}.
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
#' @details
#' New fits retain the original input dimnames in `input.dimnames`, separately
#' from custom feature IDs. Refinement and landmarks require the same dimensions
#' and reject mismatching names/order when both the fitted and supplied axis
#' are named. They never reorder automatically. Unnamed axes and older fits
#' lacking this metadata retain positional matching; names alone cannot detect
#' changed values or permutations within duplicated names.
#'
#' The `history` list records settings for each stored depth. The first entry
#' records fitting support, policy and ties. Refinement entries additionally
#' record selection mode, selected lineage/node IDs, nodes actually refined,
#' source policy, separator and refinement factor. `available = FALSE` marks
#' older levels whose history was not saved when a legacy fit is refined.
#' The top-level `n0` remains the original depth-1 threshold.
#'
#' Fitted objects store explicit feature paths and node metadata in
#' \code{nodes} and \code{node.ids}, with \code{.pure} and \code{.absorb}
#' variants. A node path contains original column indices; zero marks a
#' synthetic rare component. Node keys are local to the fitted feature mapping.
#' Metadata records parent keys, leaf-feature index, rare status and whether
#' the node was carried unchanged from the preceding stored level.
#'
#' Readable lineage IDs and labels are preserved when unambiguous. If different
#' paths produce the same string (for example, feature \code{A__B} and path
#' \code{A} then \code{B}), colliding strings receive a \code{[node ...]}
#' suffix. Use the returned IDs for explicit refinement; never split them to
#' recover features. Real features named like \code{rare.label} remain distinct
#' from synthetic rare groups. Consumers upgrade unambiguous older objects on
#' use and reject ambiguous legacy hierarchies with a rebuild message.
#'
#' @examples
#' X <- rbind(
#'   s1 = c(A = 10, B = 2, C = 1),
#'   s2 = c(A = 9, B = 3, C = 1),
#'   s3 = c(A = 1, B = 10, C = 2),
#'   s4 = c(A = 1, B = 9, C = 3),
#'   s5 = c(A = 1, B = 2, C = 10)
#' )
#' fit <- linf.csts(normalize.linf(X), n0 = 2, low.freq.policy = "pure")
#' table(fit$lineage.label)
#' fit$retained.feature.ids
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

    if (!is.numeric(n0) || length(n0) != 1L || !is.finite(n0) || n0 < 1 || n0 %% 1 != 0) {
        stop("linf.csts: n0 must be integer >= 1")
    }
    if (!is.character(rare.label) || length(rare.label) != 1L || is.na(rare.label) || !nzchar(rare.label)) {
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

                    j <- sort(kept.idx[kvals == m])
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
        input.dimnames   = dimnames(X) %||% list(NULL, NULL),
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
    out$history <- list(list(
        depth = 1L, available = TRUE, operation = "fit", n0 = as.integer(n0),
        low.freq.policy = low.freq.policy, rare.label = rare.label,
        tie.method = tie.method, backend = backend
    ))

    for (view in c("pure", "absorb")) {
        indices <- out[[paste0("depth1.feature.index.", view)]]
        paths <- lapply(indices, function(i) {
            if (!is.na(i)) as.integer(i) else if (view == "pure") 0L else integer()
        })
        out <- linf.store.node.level(out, view, paths, 1L)
    }
    out <- linf.activate.nodes(out)
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
