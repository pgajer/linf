linf.resolve.backend <- function(X, backend = c("auto", "dense", "sparse")) {
  backend <- match.arg(backend)
  if (backend == "auto") {
    if (inherits(X, "sparseMatrix")) {
      return("sparse")
    }
    return("dense")
  }
  backend
}

linf.prepare.matrix <- function(X, backend = c("auto", "dense", "sparse"), fun.name = "linf") {
  backend <- linf.resolve.backend(X, backend)

  if (backend == "dense") {
    X <- as.matrix(X)
    storage.mode(X) <- "numeric"
    return(list(X = X, backend = backend))
  }

  if (!requireNamespace("Matrix", quietly = TRUE)) {
    stop(fun.name, ": sparse backend requires the 'Matrix' package")
  }

  if (inherits(X, "sparseMatrix")) {
    X <- as(X, "dgCMatrix")
  } else {
    X <- Matrix::Matrix(X, sparse = TRUE)
    X <- as(X, "dgCMatrix")
  }

  X <- Matrix::drop0(X)
  X@x <- as.numeric(X@x)
  list(X = X, backend = backend)
}

linf.validate.matrix <- function(X, backend, fun.name) {
  if (backend == "dense") {
    if (any(!is.finite(X))) {
      stop(fun.name, ": non-finite entries found")
    }
    if (any(X < 0, na.rm = TRUE)) {
      stop(fun.name, ": negative entries found")
    }
  } else {
    if (any(!is.finite(X@x))) {
      stop(fun.name, ": non-finite entries found")
    }
    if (any(X@x < 0, na.rm = TRUE)) {
      stop(fun.name, ": negative entries found")
    }
  }

  if (nrow(X) == 0L || ncol(X) == 0L) {
    stop(fun.name, ": matrix has zero rows or columns")
  }
}

linf.row.max.sparse <- function(X) {
  Xr <- as(X, "RsparseMatrix")
  p <- Xr@p
  x <- Xr@x
  out <- numeric(nrow(Xr))

  for (ii in seq_len(nrow(Xr))) {
    start <- p[[ii]] + 1L
    end <- p[[ii + 1L]]
    if (start <= end) {
      out[[ii]] <- max(x[start:end])
    }
  }

  out
}

linf.normalize.sparse <- function(X, tol = 0) {
  Xr <- as(X, "RsparseMatrix")
  p <- Xr@p
  row.max <- linf.row.max.sparse(Xr)
  keep.rows <- which(row.max > tol)

  if (length(keep.rows)) {
    for (ii in keep.rows) {
      start <- p[[ii]] + 1L
      end <- p[[ii + 1L]]
      if (start <= end) {
        Xr@x[start:end] <- Xr@x[start:end] / row.max[[ii]]
      }
    }
  }

  Matrix::drop0(as(Xr, "CsparseMatrix"))
}

linf.choose.index <- function(candidates, tie.method, error.message) {
  if (!length(candidates)) {
    return(NA_integer_)
  }
  if (length(candidates) == 1L) {
    return(as.integer(candidates[[1L]]))
  }
  if (tie.method == "first") {
    return(as.integer(min(candidates)))
  }
  if (tie.method == "random") {
    return(as.integer(sample(candidates, 1L)))
  }
  stop(error.message)
}

linf.cells.sparse <- function(S,
                              feature.ids = NULL,
                              feature.labels = NULL,
                              tie.method = c("first", "random", "error"),
                              return.value = FALSE) {
  tie.method <- match.arg(tie.method)

  prep <- linf.prepare.matrix(S, backend = "sparse", fun.name = "linf.cells")
  X <- prep$X
  linf.validate.matrix(X, backend = "sparse", fun.name = "linf.cells")

  meta <- resolve.linf.feature.meta(X, feature.ids = feature.ids, feature.labels = feature.labels)
  id.lev <- meta$feature.ids
  lev <- meta$feature.labels

  Xr <- as(X, "RsparseMatrix")
  p <- Xr@p
  j <- Xr@j
  x <- Xr@x

  idx <- rep(NA_integer_, nrow(Xr))
  val <- numeric(nrow(Xr))

  for (ii in seq_len(nrow(Xr))) {
    start <- p[[ii]] + 1L
    end <- p[[ii + 1L]]
    if (start > end) {
      next
    }

    row.vals <- x[start:end]
    row.cols <- j[start:end] + 1L
    row.max <- max(row.vals)
    val[[ii]] <- row.max

    ties <- row.cols[row.vals == row.max]
    idx[[ii]] <- linf.choose.index(
      ties,
      tie.method = tie.method,
      error.message = "linf.cells: tie encountered and tie.method = 'error'"
    )
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

linf.absorb.sparse <- function(X,
                               to.absorb,
                               kept.idx,
                               tie.method = c("first", "random", "error"),
                               fallback.idx) {
  tie.method <- match.arg(tie.method)

  Xr <- as(X, "RsparseMatrix")
  p <- Xr@p
  j <- Xr@j
  x <- Xr@x

  kept.map <- logical(ncol(Xr))
  kept.map[kept.idx] <- TRUE
  new.idx <- integer(length(to.absorb))

  for (kk in seq_along(to.absorb)) {
    ii <- to.absorb[[kk]]
    start <- p[[ii]] + 1L
    end <- p[[ii + 1L]]

    if (start > end) {
      new.idx[[kk]] <- fallback.idx
      next
    }

    row.cols <- j[start:end] + 1L
    row.vals <- x[start:end]
    keep.mask <- kept.map[row.cols]

    if (!any(keep.mask)) {
      new.idx[[kk]] <- fallback.idx
      next
    }

    kept.cols <- row.cols[keep.mask]
    kept.vals <- row.vals[keep.mask]
    row.max <- max(kept.vals)

    if (!is.finite(row.max) || row.max <= 0) {
      new.idx[[kk]] <- fallback.idx
      next
    }

    ties <- kept.cols[kept.vals == row.max]
    new.idx[[kk]] <- linf.choose.index(
      ties,
      tie.method = tie.method,
      error.message = "linf.csts: tie during reassignment and tie.method = 'error'"
    )
  }

  new.idx
}
