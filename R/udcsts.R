#' Unordered Dominance Community State Types
#'
#' Group samples by their top-feature sets, ignoring dominance order, and then
#' apply a minimum support threshold to the pooled sets.
#'
#' @param X Nonnegative numeric samples-by-features matrix or numeric data frame.
#'   Numeric dense and sparse `Matrix` objects are also supported. Every row
#'   must contain at least one positive value.
#' @param level Integer between 1 and `ncol(X)`. Number of features in each set.
#' @param n0 Integer >= 1. Minimum pure set frequency, evaluated after pooling
#'   all dominance orders belonging to the same set.
#' @param low.freq.policy One of `"pure"`, `"absorb"`, or `"filter"`.
#'   `"pure"` marks low-support samples as rare without reassigning them;
#'   `"filter"` excludes them through the returned `keep` mask; `"absorb"`
#'   assigns them to retained sets using the rule below.
#' @param feature.ids,feature.labels Optional stable feature IDs and display
#'   labels, respectively. Each has length `ncol(X)`. Defaults and name
#'   disambiguation follow [linf.csts()].
#' @param backend One of `"auto"`, `"dense"`, or `"sparse"`. Automatic mode
#'   preserves sparse inputs. Sparse mode never densifies the whole matrix;
#'   absorption materializes only the current sample row.
#'
#' @details
#' For each row, decreasing abundance and then original column order determine
#' the ordered top `level` features. Sorting their column indices gives the
#' unordered set. Set IDs join these indices with commas, independently of
#' feature labels, so IDs are local to the fitted feature dictionary. Exact
#' ties, including zero-valued ranks, always use original column order. Rows
#' with fewer than `level` positive entries still receive a complete set;
#' `positive.features` identifies these rows. All-zero rows are rejected.
#'
#' The pure top-feature sets are formed on all input rows before thresholding.
#' For example, three `(A, B)` rows and three `(B, A)` rows form one set of size
#' six even when neither ordered group meets `n0 = 5`. At level 1, pure set
#' membership agrees with dominant-feature grouping under first-column ties.
#' Across levels the unordered partitions need not form a refinement hierarchy.
#' The result is not a `linf.csts` hierarchy and cannot be passed to its
#' refinement, policy-view, landmark, or frozen-transfer functions.
#'
#' Under `"pure"`, low-support samples have `assignment = NA` and
#' `status = "rare"`; they remain included in `keep`. Under `"filter"`,
#' they instead have `status = "filtered"` and `keep = FALSE`. Neither
#' policy alters `pure.id`, `support`, or the input matrix. An empty retained
#' family is valid for both policies. The mask selects observations, not
#' features: use `X[fit$keep, , drop = FALSE]` to obtain the filtered matrix.
#'
#' Under `"absorb"`, each low-support row goes to the retained set maximizing
#' the sum of its feature abundances. Scores are evaluated after dividing the
#' row by its positive maximum to avoid overflow. This is equivalent to
#' maximizing the unscaled sum in exact arithmetic. Ties use lexicographic
#' order of the sets' sorted numerical column-index tuples, including when
#' all candidate scores are zero. If no set survives, absorption is undefined
#' and raises an error. This symmetric set-based rule is distinct from
#' parent-local ordered dCST absorption, and can change under nonlinear
#' rank-preserving transformations even though the pure labels do not.
#'
#' Frequencies always refer to the pure sets before policy application.
#' For a set with frequency \eqn{m}, `cumm.pct` is
#' \deqn{100\sum_{B:m_B\geq m}m_B/N,}
#' where \eqn{m_B} is another pure set's frequency and \eqn{N} is the full
#' input row count. Equal-frequency rows have equal cumulative percentages.
#'
#' @return A list with:
#' \describe{
#'   \item{pure.id}{Original pure set ID for every input row.}
#'   \item{assignment}{Active set ID, or `NA` for rare/filtered rows.}
#'   \item{status}{Per-row `"retained"`, `"rare"`, `"filtered"`, or
#'     `"absorbed"`. Retained means the row's original pure set meets `n0`.}
#'   \item{keep}{Logical mask in input row order. All `TRUE` except for
#'     excluded rows under `"filter"`.}
#'   \item{support}{Table of all pure set frequencies, including small sets.}
#'   \item{frequency.table}{Three-column data frame: `ID`, `freq`, and
#'     `cumm.pct`. Sorted by decreasing frequency, then lexical set ID.}
#'   \item{retained.ids}{Set IDs meeting `n0`, in numerical tuple tie order.}
#'   \item{ordered.features,feature.sets}{Integer matrices in input row order,
#'     with `level` columns: dominance-ordered and sorted feature indices.}
#'   \item{positive.features}{Number of positive features in each input row.}
#'   \item{feature.ids,feature.labels}{Resolved column identities and labels.
#'     Index these with `feature.sets` to obtain readable set members.}
#'   \item{input.dimnames}{Original row and column names, for alignment.}
#'   \item{level,n0,low.freq.policy,matrix.backend}{Construction settings.}
#' }
#'
#' @examples
#' X <- rbind(c(A = 9, B = 5, C = 1), c(A = 8, B = 4, C = 1),
#'            c(A = 4, B = 9, C = 1), c(A = 3, B = 8, C = 1),
#'            c(A = 1, B = 3, C = 9))
#' fit <- linf.udcsts(X, level = 2, n0 = 3, low.freq.policy = "filter")
#' fit$frequency.table
#' X[fit$keep, , drop = FALSE]
#' fit$feature.labels[fit$feature.sets[1, ]]
#'
#' @seealso [linf.csts()], [linf.dominant.features()]
#' @export
linf.udcsts <- function(X, level = 2L, n0 = 50L,
                        low.freq.policy = c("pure", "absorb", "filter"),
                        feature.ids = NULL, feature.labels = NULL,
                        backend = c("auto", "dense", "sparse")) {
  low.freq.policy <- match.arg(low.freq.policy)
  numeric.input <- (is.matrix(X) && is.numeric(X) && !is.complex(X)) ||
    (is.data.frame(X) && all(vapply(X, function(x) {
      is.numeric(x) && !is.complex(x)
    }, logical(1)))) || methods::is(X, "dMatrix")
  if (!numeric.input || length(dim(X)) != 2L) {
    stop("linf.udcsts: X must be a numeric matrix or numeric data frame",
         call. = FALSE)
  }
  prep <- linf.prepare.matrix(X, backend = backend, fun.name = "linf.udcsts")
  X <- prep$X
  backend <- prep$backend
  linf.validate.matrix(X, backend, "linf.udcsts")
  if (!is.numeric(level) || is.complex(level) || length(level) != 1L || !is.finite(level) ||
      level < 1 || level > ncol(X) || level %% 1 != 0) {
    stop("linf.udcsts: level must be an integer between 1 and ncol(X)",
         call. = FALSE)
  }
  if (!is.numeric(n0) || is.complex(n0) || length(n0) != 1L || !is.finite(n0) ||
      n0 < 1 || n0 %% 1 != 0) {
    stop("linf.udcsts: n0 must be an integer >= 1", call. = FALSE)
  }
  level <- as.integer(level)
  meta <- resolve.linf.feature.meta(X, feature.ids, feature.labels)
  ranked <- linf.udcst.ranks(X, level, backend)
  paths <- ranked$paths
  sets <- t(vapply(seq_len(nrow(X)), function(i) sort(paths[i, ]),
                   integer(level)))
  # vapply simplifies a length-one result to a vector.
  sets <- matrix(sets, nrow = nrow(X), ncol = level,
                 dimnames = if (is.null(rownames(X))) NULL else list(rownames(X), NULL))
  ids <- apply(sets, 1L, paste, collapse = ",")
  support <- table(ids)
  retained.ids <- names(support)[support >= n0]
  retained.sets <- sets[match(retained.ids, ids), , drop = FALSE]
  if (length(retained.ids)) {
    ordering <- do.call(order, as.data.frame(retained.sets))
    retained.ids <- retained.ids[ordering]
    retained.sets <- retained.sets[ordering, , drop = FALSE]
  }
  is.retained <- as.vector(support[ids] >= n0)
  assignment <- ids
  status <- ifelse(is.retained, "retained", "rare")
  if (low.freq.policy != "absorb") {
    assignment[!is.retained] <- NA_character_
    if (low.freq.policy == "filter") status[!is.retained] <- "filtered"
  } else if (any(!is.retained)) {
    if (!length(retained.ids)) {
      stop("linf.udcsts: no retained state is available for absorption",
           call. = FALSE)
    }
    for (i in which(!is.retained)) {
      values <- as.numeric(X[i, ])
      values <- values / max(values)
      scores <- rowSums(matrix(values[retained.sets],
                               nrow = nrow(retained.sets), ncol = level))
      assignment[i] <- retained.ids[which.max(scores)]
    }
    status[!is.retained] <- "absorbed"
  }
  list(pure.id = ids, assignment = assignment, status = status,
       keep = if (low.freq.policy == "filter") is.retained else rep(TRUE, nrow(X)),
       support = support, frequency.table = linf.udcst.frequencies(support),
       retained.ids = retained.ids, ordered.features = paths, feature.sets = sets,
       positive.features = ranked$positive.features,
       feature.ids = meta$feature.ids, feature.labels = meta$feature.labels,
       input.dimnames = dimnames(X), level = level, n0 = n0,
       low.freq.policy = low.freq.policy, matrix.backend = backend)
}

linf.udcst.ranks <- function(X, level, backend) {
  n <- nrow(X)
  p <- ncol(X)
  paths <- matrix(NA_integer_, n, level, dimnames = if (is.null(rownames(X))) NULL else list(rownames(X), NULL))
  positive.features <- integer(n)
  if (backend == "sparse") X <- as(X, "RsparseMatrix")
  for (i in seq_len(n)) {
    if (backend == "sparse") {
      first <- X@p[i] + 1L
      last <- X@p[i + 1L]
      entries <- if (first <= last) seq.int(first, last) else integer()
      columns <- X@j[entries] + 1L
      values <- X@x[entries]
    } else {
      columns <- which(X[i, ] > 0)
      values <- X[i, columns]
    }
    positive.features[i] <- length(columns)
    if (!length(columns)) {
      stop("linf.udcsts: all-zero rows have no positive dominance context",
           call. = FALSE)
    }
    selected <- columns[order(-values, columns)]
    if (length(selected) < level) {
      selected <- c(selected, setdiff(seq_len(p), columns))
    }
    paths[i, ] <- selected[seq_len(level)]
  }
  list(paths = paths, positive.features = positive.features)
}

linf.udcst.frequencies <- function(support) {
  result <- data.frame(ID = names(support), freq = as.integer(support),
                       stringsAsFactors = FALSE)
  result <- result[order(-result$freq, result$ID), , drop = FALSE]
  tiers <- rle(result$freq)
  ends <- cumsum(tiers$lengths)
  result$cumm.pct <- rep(100 * cumsum(result$freq)[ends] / sum(result$freq),
                          tiers$lengths)
  rownames(result) <- NULL
  result
}
