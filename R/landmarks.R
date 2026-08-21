resolve.linf.lineage.labels <- function(csts,
                                        kind = c("id", "label"),
                                        view = c("active", "pure", "absorb")) {
  kind <- match.arg(kind)
  view <- match.arg(view)

  if (kind == "id") {
    active <- csts$lineage.ids
    pure <- csts$lineage.ids.pure
    absorb <- csts$lineage.ids.absorb
  } else {
    active <- csts$lineage.labels
    pure <- csts$lineage.labels.pure
    absorb <- csts$lineage.labels.absorb
  }

  switch(view,
         active = active,
         pure = pure,
         absorb = absorb)
}

resolve.linf.landmark.view <- function(csts,
                                       view = c("active", "pure", "absorb")) {
  view <- match.arg(view)
  if (view == "active") {
    return(linf.active.low.freq.view(csts$low.freq.policy))
  }
  view
}

resolve.linf.landmark.feature <- function(leaf.id,
                                          feature.ids,
                                          feature.labels,
                                          rare.label) {
  is.rare <- !is.na(leaf.id) && identical(leaf.id, rare.label)

  if (is.na(leaf.id) || is.rare) {
    return(list(
      computable = FALSE,
      is.rare = is.rare,
      index = NA_integer_,
      feature.id = NA_character_,
      feature.label = NA_character_
    ))
  }

  idx <- match(leaf.id, feature.ids)
  if (is.na(idx)) {
    idx <- match(leaf.id, feature.labels)
  }

  if (is.na(idx)) {
    return(list(
      computable = FALSE,
      is.rare = FALSE,
      index = NA_integer_,
      feature.id = NA_character_,
      feature.label = NA_character_
    ))
  }

  list(
    computable = TRUE,
    is.rare = FALSE,
    index = idx,
    feature.id = feature.ids[[idx]],
    feature.label = feature.labels[[idx]]
  )
}

choose.linf.landmark.index <- function(scores,
                                       tie.method = c("first", "random", "error"),
                                       minimize = FALSE) {
  tie.method <- match.arg(tie.method)

  best <- if (isTRUE(minimize)) min(scores) else max(scores)
  idx <- which(scores == best)

  if (length(idx) == 1L) return(idx)
  if (tie.method == "first") return(idx[1L])
  if (tie.method == "random") return(sample(idx, 1L))

  stop("linf.landmarks: tie encountered and tie.method = 'error'")
}

empty.linf.landmark.lineages <- function() {
  data.frame(
    lineage.id = character(0),
    lineage.label = character(0),
    lineage.size = integer(0),
    target.feature.id = character(0),
    target.feature.label = character(0),
    is.rare = logical(0),
    landmarks.computable = logical(0),
    stringsAsFactors = FALSE
  )
}

empty.linf.landmark.rows <- function() {
  data.frame(
    lineage.id = character(0),
    lineage.label = character(0),
    landmark.type = character(0),
    point.index = integer(0),
    point.name = character(0),
    target.feature.id = character(0),
    target.feature.label = character(0),
    observed.value = numeric(0),
    target.value = numeric(0),
    abs.deviation = numeric(0),
    stringsAsFactors = FALSE
  )
}

#' Landmark points for dCST dominance-lineages
#'
#' @description
#' Computes representative landmark points for the dominance-lineages of a
#' \code{"linf.csts"} object at a chosen depth and view.
#'
#' Landmark types are defined with respect to the leaf feature of the dCST path:
#' the last feature ID in the lineage ID path.
#' Lineages whose leaf token is \code{rare.label} are reported but skipped for
#' landmark computation because they do not correspond to a unique target
#' feature.
#'
#' @param M Numeric matrix (samples x features) used to build or refine the dCSTs.
#' @param csts A \code{"linf.csts"} object.
#' @param depth Integer. dCST depth to inspect. Defaults to the leaf depth
#'   \code{csts$depth}.
#' @param view Character. One of \code{"active"}, \code{"pure"}, or
#'   \code{"absorb"}.
#' @param landmark.types Character vector containing any of
#'   \code{"endpoint.max"}, \code{"endpoint.min"}, \code{"mean.rep"}, or
#'   \code{"median.rep"}.
#' @param tie.method Character. Tie handling for landmark selection:
#'   \code{"first"}, \code{"random"}, or \code{"error"}.
#' @param backend Character. Matrix backend to use: \code{"auto"},
#'   \code{"dense"}, or \code{"sparse"}. The default \code{"auto"} preserves
#'   sparse input and otherwise uses the dense path.
#'
#' @return A list of class \code{"linf.landmarks"} with components:
#' \itemize{
#'   \item \code{depth}, \code{view}, \code{sep}, \code{rare.label}
#'   \item \code{feature.ids}, \code{feature.labels}
#'   \item \code{lineages}: one row per dominance-lineage with computability
#'     metadata
#'   \item \code{landmarks}: one row per computed landmark point
#' }
#'
#' @export
linf.landmarks <- function(M,
                           csts,
                           depth = NULL,
                           view = c("active", "pure", "absorb"),
                           landmark.types = c("endpoint.max",
                                              "endpoint.min",
                                              "mean.rep",
                                              "median.rep"),
                           tie.method = c("first", "random", "error"),
                           backend = c("auto", "dense", "sparse")) {
  validate.linf.csts(csts)

  view <- match.arg(view)
  tie.method <- match.arg(tie.method)
  landmark.types <- match.arg(landmark.types, several.ok = TRUE)
  if (missing(backend)) {
    backend <- csts$matrix.backend
  }

  prep <- linf.prepare.matrix(M, backend = backend, fun.name = "linf.landmarks")
  X <- prep$X
  backend <- prep$backend
  linf.validate.matrix(X, backend = backend, fun.name = "linf.landmarks")
  if (nrow(X) != length(csts$lineage.label)) {
    stop("linf.landmarks: nrow(M) must match the number of samples in csts")
  }

  max.depth <- csts$depth
  if (is.null(depth)) depth <- max.depth

  if (!is.numeric(depth) || length(depth) != 1L || depth < 1L || depth %% 1 != 0) {
    stop("linf.landmarks: depth must be an integer >= 1")
  }
  if (depth > max.depth) {
    stop("linf.landmarks: requested depth exceeds csts$depth")
  }

  resolved.view <- resolve.linf.landmark.view(csts, view)
  id.levels <- resolve.linf.lineage.labels(csts, kind = "id", view = resolved.view)
  label.levels <- resolve.linf.lineage.labels(csts, kind = "label", view = resolved.view)

  if (length(id.levels) < depth || length(label.levels) < depth) {
    stop("linf.landmarks: requested view does not contain the requested depth")
  }

  lineage.ids <- as.character(id.levels[[depth]])
  lineage.labels <- as.character(label.levels[[depth]])

  if (length(lineage.ids) != nrow(X) || length(lineage.labels) != nrow(X)) {
    stop("linf.landmarks: dCST levels at the requested depth must align with rows of M")
  }

  meta <- resolve.linf.feature.meta(
    X,
    feature.ids = csts$feature.ids,
    feature.labels = csts$feature.labels
  )

  sep <- csts$sep %||% "__"
  rare.label <- csts$rare.label
  row.ids <- rownames(X)
  if (is.null(row.ids)) row.ids <- rep(NA_character_, nrow(X))

  lineages <- empty.linf.landmark.lineages()
  landmarks <- empty.linf.landmark.rows()
  unique.lineages <- unique(lineage.ids[!is.na(lineage.ids)])

  for (lineage in unique.lineages) {
    members <- which(lineage.ids == lineage)
    lineage.label <- lineage.labels[members[1L]]
    parts <- strsplit(lineage, sep, fixed = TRUE)[[1L]]
    leaf.id <- parts[[length(parts)]]
    target <- resolve.linf.landmark.feature(
      leaf.id,
      feature.ids = meta$feature.ids,
      feature.labels = meta$feature.labels,
      rare.label = rare.label
    )

    lineages <- rbind(lineages, data.frame(
      lineage.id = lineage,
      lineage.label = lineage.label,
      lineage.size = length(members),
      target.feature.id = target$feature.id,
      target.feature.label = target$feature.label,
      is.rare = target$is.rare,
      landmarks.computable = target$computable,
      stringsAsFactors = FALSE
    ))

    if (!isTRUE(target$computable)) next

    vals <- as.numeric(X[members, target$index, drop = FALSE])

    for (landmark.type in landmark.types) {
      if (landmark.type == "endpoint.max") {
        pick <- choose.linf.landmark.index(vals, tie.method = tie.method, minimize = FALSE)
        observed <- vals[[pick]]
        target.value <- observed
      } else if (landmark.type == "endpoint.min") {
        pick <- choose.linf.landmark.index(vals, tie.method = tie.method, minimize = TRUE)
        observed <- vals[[pick]]
        target.value <- observed
      } else if (landmark.type == "mean.rep") {
        target.value <- mean(vals)
        dev <- abs(vals - target.value)
        pick <- choose.linf.landmark.index(dev, tie.method = tie.method, minimize = TRUE)
        observed <- vals[[pick]]
      } else if (landmark.type == "median.rep") {
        target.value <- stats::median(vals)
        dev <- abs(vals - target.value)
        pick <- choose.linf.landmark.index(dev, tie.method = tie.method, minimize = TRUE)
        observed <- vals[[pick]]
      } else {
        stop("linf.landmarks: unsupported landmark type")
      }

      point.index <- members[[pick]]
      landmarks <- rbind(landmarks, data.frame(
        lineage.id = lineage,
        lineage.label = lineage.label,
        landmark.type = landmark.type,
        point.index = point.index,
        point.name = row.ids[[point.index]],
        target.feature.id = target$feature.id,
        target.feature.label = target$feature.label,
        observed.value = observed,
        target.value = target.value,
        abs.deviation = abs(observed - target.value),
        stringsAsFactors = FALSE
      ))
    }
  }

  structure(
    list(
      depth = as.integer(depth),
      view = resolved.view,
      sep = sep,
      rare.label = rare.label,
      feature.ids = meta$feature.ids,
      feature.labels = meta$feature.labels,
      lineages = lineages,
      landmarks = landmarks
    ),
    class = "linf.landmarks"
  )
}
