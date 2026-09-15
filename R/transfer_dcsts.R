#' Transfer samples into a frozen dCST hierarchy
#'
#' @description
#' Assigns rows of a new sample-by-feature matrix into the realized lineage sets
#' of a fitted \code{"linf.csts"} hierarchy. The hierarchy is not refit: each
#' sample is walked through the frozen tree by choosing, at each depth, the
#' realized child lineage whose newly added feature has the largest abundance in
#' that sample. Scoring uses explicit node feature indices, independently of
#' display labels and separator-containing feature IDs.
#'
#' @param X Nonnegative sample-by-feature matrix to transfer. Columns may be in
#'   any order and may include features not present in \code{csts}; missing
#'   retained features are treated as zero.
#' @param csts A fitted \code{"linf.csts"} object produced by
#'   \code{\link{linf.csts}} and optionally refined by
#'   repeated calls to \code{\link{refine.linf.csts}}.
#' @param depth Integer vector of requested depths. Defaults to all fitted
#'   depths in \code{csts}.
#' @param view Which fitted hierarchy view to transfer into. \code{"absorb"}
#'   uses \code{lineage.labels.absorb}; \code{"pure"} uses
#'   \code{lineage.labels.pure}; \code{"active"} uses \code{lineage.labels}.
#' @param match.by Whether columns in \code{X} are aligned to
#'   \code{csts$feature.ids} or \code{csts$feature.labels}.
#' @param feature.ids Optional feature identifiers for columns of \code{X}.
#'   Defaults to \code{colnames(X)}. When used for matching, keys must be
#'   unique, nonmissing and nonempty; ambiguous query keys are rejected.
#' @param feature.labels Optional feature labels for columns of \code{X}.
#'   Defaults to \code{feature.ids}. The same key requirements apply when
#'   matching by labels; repeated display labels are allowed when matching by IDs.
#' @param tie.method How to resolve ties among frozen child lineages with equal
#'   sample abundance. \code{"support"} chooses the tied lineage with largest
#'   reference support and then lexical order; \code{"first"} uses frozen child
#'   order; \code{"random"} samples one tied lineage; \code{"error"} stops.
#' @param carry.forward.terminal.depths Logical, retained for call compatibility.
#'   Under either setting, terminal lineages are assigned only where they occur
#'   in the fitted hierarchy. No out-of-hierarchy fallback is performed.
#' @param sep Separator recorded with the result for call compatibility. Node
#'   identity and feature lookup use the fitted metadata, not this string.
#'   Defaults to \code{csts$sep} or
#'   \code{"__"}.
#' @param backend Matrix backend; passed to the internal matrix preparer.
#'
#' @return A list of class \code{"linf.dcst.transfer"} with components:
#'   \itemize{
#'     \item \code{assignment}: character matrix of transferred labels for the
#'       requested depths.
#'     \item \code{all.depths}: character matrix for all fitted depths.
#'     \item \code{assignment.ids}, \code{all.depths.ids}: corresponding
#'       lineage-ID matrices whose identifiers do not encode display labels. IDs belong
#'       to the fitted feature mapping; they are not universal across fits.
#'     \item \code{feature.match}: one row per fitted feature, with ID, label,
#'       query column index (NA if absent) and a logical matched flag.
#'     \item \code{diagnostics}: one row per query, with sample index/name,
#'       deepest assigned depth, stopping depth/reason, counts of present and
#'       absent candidate features at that stop, and number of depths with ties.
#'       Diagnostics describe the full walk, regardless of requested depths.
#'     \item \code{depth}: requested depth vector.
#'     \item \code{view}, \code{match.by}, \code{tie.method}: settings used.
#'   }
#'   At any depth with no realized candidate or no positive abundance for its
#'   candidate features, that depth and all subsequent depths are \code{NA}.
#'   The original assignment matrices retain display labels regardless of
#'   \code{match.by}; use their `.ids` counterparts for joins.
#'
#' @details
#' Diagnostic reasons are `complete` (every fitted depth assigned),
#' `no_candidates` (no realized children), `synthetic_only` (only unscorable
#' rare categories), `missing_candidate_features` (all real candidate features
#' absent from the query), or `no_positive_candidate_values` (some real candidate
#' features are present but none has positive abundance). A partially missing
#' candidate set can have the last reason; inspect `missing.candidates` too.
#' Completed walks have NA stopping depth and candidate counts. A terminal
#' lineage carried through later stored levels counts as an assignment at each
#' such level; `assigned.depth` is not the number of distinct feature choices.
#' ID output does not change the existing tie rule: lexical display-label order
#' still resolves equal-support ties. Relabeling and refitting can therefore
#' change tied assignments; inspect `n.tied.depths` when comparing results.
#'
#' @examples
#' X <- rbind(
#'   s1 = c(A = 10, B = 2, C = 1),
#'   s2 = c(A = 9, B = 3, C = 1),
#'   s3 = c(A = 1, B = 10, C = 2),
#'   s4 = c(A = 1, B = 9, C = 3)
#' )
#' M <- normalize.linf(X)
#' fit <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")
#' transfer.dcsts(X, fit)$assignment
#'
#' @export
transfer.dcsts <- function(X,
                           csts,
                           depth = NULL,
                           view = c("absorb", "active", "pure"),
                           match.by = c("feature.ids", "feature.labels"),
                           feature.ids = NULL,
                           feature.labels = NULL,
                           tie.method = c("support", "first", "random", "error"),
                           carry.forward.terminal.depths = TRUE,
                           sep = NULL,
                           backend = c("auto", "dense", "sparse")) {
  validate.linf.csts(csts)
  csts <- linf.ensure.nodes(csts)
  view <- match.arg(view)
  match.by <- match.arg(match.by)
  tie.method <- match.arg(tie.method)
  if (is.null(sep)) sep <- csts$sep %||% "__"
  if (!is.character(sep) || length(sep) != 1L || is.na(sep) || !nzchar(sep)) {
    stop("transfer.dcsts: sep must be a non-empty character scalar")
  }
  if (!is.logical(carry.forward.terminal.depths) ||
      length(carry.forward.terminal.depths) != 1L ||
      is.na(carry.forward.terminal.depths)) {
    stop("transfer.dcsts: carry.forward.terminal.depths must be TRUE or FALSE")
  }

  prep <- linf.prepare.matrix(X, backend = backend, fun.name = "transfer.dcsts")
  X <- prep$X
  backend <- prep$backend
  linf.validate.matrix(X, backend = backend, fun.name = "transfer.dcsts")

  feature.ids <- feature.ids %||% colnames(X)
  if (is.null(feature.ids)) feature.ids <- paste0("V", seq_len(ncol(X)))
  feature.labels <- feature.labels %||% feature.ids
  if (length(feature.ids) != ncol(X) || length(feature.labels) != ncol(X)) {
    stop("transfer.dcsts: feature.ids and feature.labels must have length ncol(X)")
  }
  # Only the selected matching vector must be unique; display labels may repeat
  # when stable IDs are the selected keys.
  if (match.by == "feature.ids") {
    feature.ids <- linf.validate.query.keys(feature.ids, "feature.ids", ncol(X))
  } else {
    feature.labels <- linf.validate.query.keys(feature.labels, "feature.labels", ncol(X))
  }

  resolved.view <- resolve.linf.landmark.view(csts, view)
  levels <- csts[[paste0("node.ids.", resolved.view)]]
  nodes <- csts[[paste0("nodes.", resolved.view)]]
  max.depth <- length(levels)
  if (!max.depth) stop("transfer.dcsts: csts does not contain fitted levels")
  if (is.null(depth)) depth <- seq_len(max.depth)
  if (!is.numeric(depth) || any(!is.finite(depth)) || any(depth < 1) ||
      any(depth %% 1 != 0) || any(depth > max.depth)) {
    stop("transfer.dcsts: depth must contain fitted depths between 1 and csts$depth")
  }
  depth <- as.integer(depth)

  ref.features <- switch(
    match.by,
    feature.ids = csts$feature.ids,
    feature.labels = csts$feature.labels
  )
  ref.features <- as.character(ref.features)

  query.features <- switch(
    match.by,
    feature.ids = as.character(feature.ids),
    feature.labels = as.character(feature.labels)
  )
  X.aligned <- linf.align.transfer.matrix(X, query.features, ref.features, backend = backend)
  # Alignment and node scoring use feature positions; labels are presentation.
  colnames(X.aligned) <- csts$feature.labels

  tree <- linf.dcst.transfer.tree(levels, nodes, max.depth = max.depth, sep = sep)
  query.index <- match(ref.features, query.features)
  matched <- !is.na(query.index)
  walks <- lapply(seq_len(nrow(X.aligned)), function(i) {
    linf.transfer.one.sample(
      as.numeric(X.aligned[i, , drop = TRUE]), tree = tree,
      max.depth = max.depth, tie.method = tie.method, sep = sep,
      feature.matched = matched
    )
  })
  node.assignments <- do.call(rbind, lapply(walks, `[[`, "assignment"))
  rownames(node.assignments) <- rownames(X)
  colnames(node.assignments) <- paste0("depth", seq_len(max.depth))
  assignments <- assignment.ids <- node.assignments
  for (d in seq_len(max.depth)) {
    index <- match(node.assignments[, d], nodes[[d]]$node.id)
    assignments[, d] <- nodes[[d]]$lineage.label[index]
    assignment.ids[, d] <- nodes[[d]]$lineage.id[index]
  }
  diagnostics <- do.call(rbind, lapply(walks, `[[`, "diagnostics"))
  diagnostics <- data.frame(sample.index = seq_len(nrow(X)),
                            sample.id = rownames(X) %||% rep(NA_character_, nrow(X)),
                            diagnostics, row.names = NULL)
  feature.match <- data.frame(
    feature.id = csts$feature.ids, feature.label = csts$feature.labels,
    query.index = query.index, matched = matched, stringsAsFactors = FALSE
  )

  out <- list(
    assignment = assignments[, depth, drop = FALSE],
    all.depths = assignments,
    assignment.ids = assignment.ids[, depth, drop = FALSE],
    all.depths.ids = assignment.ids,
    diagnostics = diagnostics,
    feature.match = feature.match,
    depth = depth,
    view = view,
    match.by = match.by,
    tie.method = tie.method,
    carry.forward.terminal.depths = carry.forward.terminal.depths,
    sep = sep
  )
  class(out) <- "linf.dcst.transfer"
  out
}

linf.align.transfer.matrix <- function(X, query.features, ref.features, backend) {
  idx <- match(ref.features, query.features)
  if (backend == "sparse") {
    out <- Matrix::Matrix(0, nrow = nrow(X), ncol = length(ref.features), sparse = TRUE)
    hit <- !is.na(idx)
    if (any(hit)) out[, hit] <- X[, idx[hit], drop = FALSE]
    rownames(out) <- rownames(X)
    return(out)
  }
  out <- matrix(0, nrow = nrow(X), ncol = length(ref.features))
  hit <- !is.na(idx)
  if (any(hit)) out[, hit] <- X[, idx[hit], drop = FALSE]
  rownames(out) <- rownames(X)
  out
}

linf.dcst.transfer.tree <- function(levels, nodes, max.depth, sep) {
  tree <- vector("list", max.depth)
  support <- vector("list", max.depth)

  level1 <- as.character(levels[[1L]])
  level1 <- level1[!is.na(level1) & nzchar(level1)]
  root.counts <- table(level1)
  root.labels <- nodes[[1L]]$lineage.label[match(names(root.counts), nodes[[1L]]$node.id)]
  root.counts <- root.counts[order(-as.numeric(root.counts), root.labels)]
  tree[[1L]] <- list("__ROOT__" = names(root.counts))
  support[[1L]] <- list("__ROOT__" = as.numeric(root.counts))
  names(support[[1L]][["__ROOT__"]]) <- names(root.counts)

  if (max.depth >= 2L) {
    for (d in 2:max.depth) {
      parent <- as.character(levels[[d - 1L]])
      child <- as.character(levels[[d]])
      keep <- !is.na(parent) & nzchar(parent) & !is.na(child) & nzchar(child)
      parent <- parent[keep]
      child <- child[keep]
      pairs <- unique(data.frame(parent = parent, child = child, stringsAsFactors = FALSE))
      counts <- as.data.frame(table(parent, child), stringsAsFactors = FALSE)
      counts <- counts[counts$Freq > 0, , drop = FALSE]

      tree[[d]] <- list()
      support[[d]] <- list()
      if (!nrow(pairs)) next
      parent.groups <- split(pairs$child, pairs$parent)
      support.groups <- split(counts, counts$parent)
      for (key in sort(names(parent.groups))) {
        children <- unique(parent.groups[[key]])
        tree[[d]][[key]] <- children
        counts.df <- support.groups[[key]]
        counts.vec <- stats::setNames(counts.df$Freq, as.character(counts.df$child))
        support[[d]][[key]] <- counts.vec
      }
    }
  }

  list(children = tree, support = support, nodes = nodes, sep = sep)
}

linf.transfer.one.sample <- function(sample.values,
                                     tree,
                                     max.depth,
                                     tie.method,
                                     sep,
                                     feature.matched) {
  out <- rep(NA_character_, max.depth)
  parent.key <- "__ROOT__"
  reason <- "complete"
  stop.depth <- NA_integer_
  matched.at.stop <- missing.at.stop <- NA_integer_
  tie.depths <- integer()
  child.score <- function(key, depth) {
    nodes <- tree$nodes[[depth]]
    node <- match(key, nodes$node.id)
    if (nodes$is.rare[[node]]) return(0)
    as.numeric(sample.values[[nodes$feature.index[[node]]]])
  }

  for (d in seq_len(max.depth)) {
    candidates <- tree$children[[d]][[parent.key]]
    if (is.null(candidates) || !length(candidates)) {
      reason <- "no_candidates"
      stop.depth <- d
      matched.at.stop <- missing.at.stop <- 0L
      break
    }

    vals <- vapply(candidates, child.score, numeric(1L), depth = d)
    best <- max(vals)
    if (!is.finite(best) || best <= 0) {
      indices <- tree$nodes[[d]]$feature.index[match(candidates, tree$nodes[[d]]$node.id)]
      indices <- unique(indices[!is.na(indices)])
      matched.at.stop <- sum(feature.matched[indices])
      missing.at.stop <- sum(!feature.matched[indices])
      reason <- if (!length(indices)) "synthetic_only" else if (!matched.at.stop) {
        "missing_candidate_features"
      } else "no_positive_candidate_values"
      stop.depth <- d
      break
    }

    top <- candidates[vals == best]
    if (length(top) > 1L) {
      tie.depths <- c(tie.depths, d)
      if (tie.method == "support") {
        supp <- tree$support[[d]][[parent.key]][top]
        supp[is.na(supp)] <- 0
        top <- top[supp == max(supp)]
        if (length(top) > 1L) {
          labels <- tree$nodes[[d]]$lineage.label[match(top, tree$nodes[[d]]$node.id)]
          top <- top[order(labels)]
        }
      } else if (tie.method == "first") {
        top <- top[[1L]]
      } else if (tie.method == "random") {
        top <- sample(top, 1L)
      } else {
        stop("transfer.dcsts: tie during frozen transfer and tie.method = 'error'")
      }
    }

    chosen <- top[[1L]]
    out[[d]] <- chosen
    parent.key <- chosen
  }

  list(assignment = out, diagnostics = data.frame(
    assigned.depth = sum(!is.na(out)), stop.depth = stop.depth,
    reason = reason, matched.candidates = matched.at.stop,
    missing.candidates = missing.at.stop, n.tied.depths = length(tie.depths),
    stringsAsFactors = FALSE
  ))
}
