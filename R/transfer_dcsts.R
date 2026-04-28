#' Transfer samples into a frozen dCST hierarchy
#'
#' @description
#' Assigns rows of a new sample-by-feature matrix into the realized lineage sets
#' of a fitted \code{"linf.csts"} hierarchy. The hierarchy is not refit: each
#' sample is walked through the frozen tree by choosing, at each depth, the
#' realized child lineage whose newly added feature has the largest abundance in
#' that sample.
#'
#' @param X Nonnegative sample-by-feature matrix to transfer. Columns may be in
#'   any order and may include features not present in \code{csts}; missing
#'   retained features are treated as zero.
#' @param csts A fitted \code{"linf.csts"} object produced by
#'   \code{\link{linf.csts}} and optionally refined by
#'   \code{\link{refine.linf.csts}} or \code{\link{refine.linf.csts.iter}}.
#' @param depth Integer vector of requested depths. Defaults to all fitted
#'   depths in \code{csts}.
#' @param view Which fitted hierarchy view to transfer into. \code{"absorb"}
#'   uses \code{cst.levels.absorb}; \code{"rare"} uses
#'   \code{cst.levels.rare}; \code{"active"} uses \code{cst.levels}.
#' @param match.by Whether columns in \code{X} are aligned to
#'   \code{csts$feature.ids} or \code{csts$feature.labels}.
#' @param feature.ids Optional feature identifiers for columns of \code{X}.
#'   Defaults to \code{colnames(X)}.
#' @param feature.labels Optional feature labels for columns of \code{X}.
#'   Defaults to \code{feature.ids}.
#' @param tie.method How to resolve ties among frozen child lineages with equal
#'   sample abundance. \code{"support"} chooses the tied lineage with largest
#'   reference support and then lexical order; \code{"first"} uses frozen child
#'   order; \code{"random"} samples one tied lineage; \code{"error"} stops.
#' @param carry.forward.terminal.depths Logical. If \code{TRUE}, an unmatched
#'   depth-1 sample may fall back to its dominant retained feature. Terminal
#'   lineages at later depths are assigned only when they appear as realized
#'   stable lineage sets in the fitted hierarchy.
#' @param sep Separator used in lineage labels. Defaults to \code{csts$sep} or
#'   \code{"__"}.
#' @param backend Matrix backend; passed to the internal matrix preparer.
#'
#' @return A list of class \code{"linf.dcst.transfer"} with components:
#'   \itemize{
#'     \item \code{assignment}: character matrix of transferred labels for the
#'       requested depths.
#'     \item \code{all.depths}: character matrix for all fitted depths.
#'     \item \code{depth}: requested depth vector.
#'     \item \code{view}, \code{match.by}, \code{tie.method}: settings used.
#'   }
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
                           view = c("absorb", "active", "rare"),
                           match.by = c("feature.ids", "feature.labels"),
                           feature.ids = NULL,
                           feature.labels = NULL,
                           tie.method = c("support", "first", "random", "error"),
                           carry.forward.terminal.depths = TRUE,
                           sep = NULL,
                           backend = c("auto", "dense", "sparse")) {
  validate.linf.csts(csts)
  view <- match.arg(view)
  match.by <- match.arg(match.by)
  tie.method <- match.arg(tie.method)
  if (is.null(sep)) sep <- csts$sep %||% "__"
  if (!is.character(sep) || length(sep) != 1L || !nzchar(sep)) {
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
  if (is.null(feature.ids)) feature.ids <- paste0("feature", seq_len(ncol(X)))
  feature.labels <- feature.labels %||% feature.ids
  if (length(feature.ids) != ncol(X)) {
    stop("transfer.dcsts: feature.ids must have length ncol(X)")
  }
  if (length(feature.labels) != ncol(X)) {
    stop("transfer.dcsts: feature.labels must have length ncol(X)")
  }

  levels <- resolve.linf.cst.levels(csts, kind = "label", view = view)
  max.depth <- length(levels)
  if (!max.depth) stop("transfer.dcsts: csts does not contain fitted levels")
  if (is.null(depth)) depth <- seq_len(max.depth)
  if (!is.numeric(depth) || any(!is.finite(depth)) || any(depth < 1) ||
      any(depth %% 1 != 0) || any(depth > max.depth)) {
    stop("transfer.dcsts: depth must contain fitted depths between 1 and csts$cst.depth")
  }
  depth <- as.integer(depth)

  ref.features <- switch(
    match.by,
    feature.ids = csts$feature.ids,
    feature.labels = csts$feature.labels
  )
  if (is.null(ref.features)) {
    ref.features <- switch(
      match.by,
      feature.ids = csts$feature.labels,
      feature.labels = csts$feature.ids
    )
  }
  if (is.null(ref.features)) {
    stop("transfer.dcsts: csts does not contain feature metadata for matching")
  }
  ref.features <- as.character(ref.features)

  query.features <- switch(
    match.by,
    feature.ids = as.character(feature.ids),
    feature.labels = as.character(feature.labels)
  )
  X.aligned <- linf.align.transfer.matrix(X, query.features, ref.features, backend = backend)
  colnames(X.aligned) <- ref.features

  tree <- linf.dcst.transfer.tree(levels, max.depth = max.depth, sep = sep)
  assigned.raw <- vapply(
    seq_len(nrow(X.aligned)),
    function(i) {
      vals <- X.aligned[i, , drop = TRUE]
      names(vals) <- ref.features
      linf.transfer.one.sample(
        vals,
        tree = tree,
        max.depth = max.depth,
        tie.method = tie.method,
        carry.forward.terminal.depths = carry.forward.terminal.depths,
        sep = sep
      )
    },
    character(max.depth)
  )
  assignments <- if (is.null(dim(assigned.raw))) {
    matrix(assigned.raw, ncol = 1L)
  } else {
    t(assigned.raw)
  }
  rownames(assignments) <- rownames(X)
  colnames(assignments) <- paste0("depth", seq_len(max.depth))

  out <- list(
    assignment = assignments[, depth, drop = FALSE],
    all.depths = assignments,
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

linf.dcst.transfer.tree <- function(levels, max.depth, sep) {
  tree <- vector("list", max.depth)
  support <- vector("list", max.depth)

  level1 <- as.character(levels[[1L]])
  level1 <- level1[!is.na(level1) & nzchar(level1)]
  root.counts <- sort(table(level1), decreasing = TRUE)
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

  list(children = tree, support = support, sep = sep)
}

linf.transfer.one.sample <- function(sample.values,
                                     tree,
                                     max.depth,
                                     tie.method,
                                     carry.forward.terminal.depths,
                                     sep) {
  out <- rep(NA_character_, max.depth)
  parent.key <- "__ROOT__"

  fallback.to.dominant <- function() {
    vals <- sample.values
    vals[is.na(vals)] <- 0
    best <- max(vals)
    if (!is.finite(best) || best <= 0) return(NA_character_)
    top <- names(vals)[vals == best]
    sort(top)[[1L]]
  }

  child.score <- function(label) {
    component <- linf.lineage.last.component(label, sep = sep)
    if (!nzchar(component)) return(0)
    val <- unname(sample.values[component])
    if (!length(val) || is.na(val)) return(0)
    as.numeric(val[[1L]])
  }

  for (d in seq_len(max.depth)) {
    candidates <- tree$children[[d]][[parent.key]]
    if (is.null(candidates) || !length(candidates)) {
      if (d == 1L && carry.forward.terminal.depths) {
        fallback <- fallback.to.dominant()
        if (!is.na(fallback)) out[[d]] <- fallback
      }
      break
    }

    vals <- vapply(candidates, child.score, numeric(1L))
    best <- max(vals)
    if (!is.finite(best) || best <= 0) {
      if (d == 1L && carry.forward.terminal.depths) {
        fallback <- fallback.to.dominant()
        if (!is.na(fallback)) out[[d]] <- fallback
      }
      break
    }

    top <- candidates[vals == best]
    if (length(top) > 1L) {
      if (tie.method == "support") {
        supp <- tree$support[[d]][[parent.key]][top]
        supp[is.na(supp)] <- 0
        top <- top[supp == max(supp)]
        if (length(top) > 1L) top <- sort(top)
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

  out
}

linf.lineage.last.component <- function(label, sep) {
  if (is.na(label) || !nzchar(label)) return("")
  parts <- strsplit(label, sep, fixed = TRUE)[[1L]]
  utils::tail(parts, 1L)[[1L]]
}
