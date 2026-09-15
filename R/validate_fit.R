validate.linf.csts <- function(obj) {
  if (!inherits(obj, "linf.csts") || !is.list(obj)) {
    stop("dCST object must inherit from class 'linf.csts'", call. = FALSE)
  }

  required <- c(
    "depth1.feature.index", "depth1.feature.id", "depth1.feature.label",
    "depth1.feature.index.pure", "depth1.feature.id.pure",
    "depth1.feature.label.pure", "depth1.feature.index.absorb",
    "depth1.feature.id.absorb", "depth1.feature.label.absorb",
    "lineage.id", "lineage.label", "lineage.id.pure",
    "lineage.label.pure", "lineage.id.absorb", "lineage.label.absorb",
    "lineage.ids", "lineage.labels", "lineage.ids.pure",
    "lineage.labels.pure", "lineage.ids.absorb", "lineage.labels.absorb",
    "feature.ids", "feature.labels", "matrix.backend", "depth", "n0",
    "low.freq.policy", "rare.label"
  )
  missing.fields <- setdiff(required, names(obj))
  if (length(missing.fields)) {
    stop(
      paste0(
        "incompatible dCST object; rebuild it with linf >= 0.2.0. Missing: ",
        paste(missing.fields, collapse = ", ")
      ),
      call. = FALSE
    )
  }

  if (!is.numeric(obj$depth) || length(obj$depth) != 1L ||
      is.na(obj$depth) || obj$depth < 1L || obj$depth %% 1 != 0) {
    stop("dCST object has an invalid depth", call. = FALSE)
  }
  depth <- as.integer(obj$depth)
  n <- length(obj$lineage.label)

  hierarchy.fields <- c(
    "lineage.ids", "lineage.labels", "lineage.ids.pure",
    "lineage.labels.pure", "lineage.ids.absorb", "lineage.labels.absorb"
  )
  for (field in hierarchy.fields) {
    hierarchy <- obj[[field]]
    if (!is.list(hierarchy) || length(hierarchy) != depth ||
        any(lengths(hierarchy) != n)) {
      stop(
        sprintf("dCST object has an invalid %s hierarchy", field),
        call. = FALSE
      )
    }
  }

  sample.fields <- c(
    "depth1.feature.index", "depth1.feature.id", "depth1.feature.label",
    "depth1.feature.index.pure", "depth1.feature.id.pure",
    "depth1.feature.label.pure", "depth1.feature.index.absorb",
    "depth1.feature.id.absorb", "depth1.feature.label.absorb",
    "lineage.id", "lineage.label", "lineage.id.pure",
    "lineage.label.pure", "lineage.id.absorb", "lineage.label.absorb"
  )
  if (any(vapply(obj[sample.fields], length, integer(1)) != n)) {
    stop("dCST object contains sample assignments of inconsistent length", call. = FALSE)
  }

  if (length(obj$feature.ids) != length(obj$feature.labels)) {
    stop("dCST object has inconsistent feature metadata", call. = FALSE)
  }
  if (!obj$low.freq.policy %in% c("pure", "absorb")) {
    stop("dCST object has an invalid low.freq.policy", call. = FALSE)
  }

  invisible(TRUE)
}
