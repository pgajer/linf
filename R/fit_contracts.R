# Dimnames describe the fitted input, independently of custom feature IDs.
# Unnamed axes and legacy fits retain their documented positional contract.
linf.validate.fit.matrix <- function(X, csts, fun.name) {
  expected <- c(length(csts$lineage.id), length(csts$feature.ids))
  axes <- c("row", "column")
  for (axis in seq_len(2L)) {
    if (dim(X)[[axis]] != expected[[axis]]) {
      stop(fun.name, ": ", axes[[axis]], " count must match the fitted matrix (expected ",
           expected[[axis]], "); use all fitted samples and features in fitted order", call. = FALSE)
    }
    stored <- csts$input.dimnames[[axis]]
    supplied <- dimnames(X)[[axis]]
    if (!is.null(stored) && !is.null(supplied) && !identical(stored, supplied)) {
      stop(fun.name, ": ", axes[[axis]], " names/order differ from the fitted matrix; ",
           "align this axis to csts$input.dimnames[[", axis, "]] before calling", call. = FALSE)
    }
  }
  invisible(TRUE)
}

linf.validate.query.keys <- function(keys, name, n) {
  if (length(keys) != n) stop("transfer.dcsts: ", name, " must have length ncol(X)", call. = FALSE)
  keys <- as.character(keys)
  if (anyNA(keys) || any(!nzchar(keys))) {
    stop("transfer.dcsts: ", name, " must not contain missing or empty keys", call. = FALSE)
  }
  if (anyDuplicated(keys)) {
    stop("transfer.dcsts: ", name, " must be unique; resolve duplicate query keys before transfer", call. = FALSE)
  }
  keys
}

linf.fit.history <- function(csts) {
  if (!is.null(csts$history)) return(csts$history)
  lapply(seq_len(csts$depth), function(depth) list(depth = depth, available = FALSE))
}
