resolve.linf.cst.levels <- function(csts,
                                    kind = c("id", "label"),
                                    view = c("active", "rare", "absorb")) {
  kind <- match.arg(kind)
  view <- match.arg(view)

  if (kind == "id") {
    active <- csts$cst.id.levels
    if (is.null(active)) active <- list(level1 = csts$cell.id %||% csts$cell.label)

    rare <- csts$cst.id.levels.rare
    if (is.null(rare) && !is.null(csts$cell.id.rare)) rare <- list(level1 = csts$cell.id.rare)

    absorb <- csts$cst.id.levels.absorb
    if (is.null(absorb) && !is.null(csts$cell.id.absorb)) absorb <- list(level1 = csts$cell.id.absorb)
  } else {
    active <- csts$cst.levels
    if (is.null(active)) active <- list(level1 = csts$cell.label)

    rare <- csts$cst.levels.rare
    if (is.null(rare) && !is.null(csts$cell.label.rare)) rare <- list(level1 = csts$cell.label.rare)

    absorb <- csts$cst.levels.absorb
    if (is.null(absorb) && !is.null(csts$cell.label.absorb)) absorb <- list(level1 = csts$cell.label.absorb)
  }

  rare <- rare %||% active
  absorb <- absorb %||% active

  switch(view,
         active = active,
         rare = rare,
         absorb = absorb)
}

resolve.linf.landmark.view <- function(csts,
                                       view = c("active", "rare", "absorb")) {
  view <- match.arg(view)
  if (view == "active") {
    return(csts$low.freq.policy %||% "active")
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

empty.linf.landmark.cells <- function() {
  data.frame(
    cell.id = character(0),
    cell.label = character(0),
    cell.size = integer(0),
    target.feature.id = character(0),
    target.feature.label = character(0),
    is.rare = logical(0),
    landmarks.computable = logical(0),
    stringsAsFactors = FALSE
  )
}

empty.linf.landmark.rows <- function() {
  data.frame(
    cell.id = character(0),
    cell.label = character(0),
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

#' Landmark points for L-infinity CST cells
#'
#' @description
#' Computes representative landmark points for the cells of a
#' \code{"linf.csts"} object at a chosen depth and view.
#'
#' In the first implementation pass, landmark types are defined with respect to
#' the leaf feature of the CST path: the last feature ID in the cell ID path.
#' Cells whose leaf token is \code{rare.label} are reported but skipped for
#' landmark computation because they do not correspond to a unique target
#' feature.
#'
#' @param M Numeric matrix (samples x features) used to build or refine the CSTs.
#' @param csts A \code{"linf.csts"} object.
#' @param depth Integer. CST depth to inspect. Defaults to the leaf depth
#'   \code{csts$cst.depth}.
#' @param view Character. One of \code{"active"}, \code{"rare"}, or
#'   \code{"absorb"}.
#' @param landmark.types Character vector containing any of
#'   \code{"endpoint.max"}, \code{"endpoint.min"}, \code{"mean.rep"}, or
#'   \code{"median.rep"}.
#' @param tie.method Character. Tie handling for landmark selection:
#'   \code{"first"}, \code{"random"}, or \code{"error"}.
#'
#' @return A list of class \code{"linf.landmarks"} with components:
#' \itemize{
#'   \item \code{depth}, \code{view}, \code{sep}, \code{rare.label}
#'   \item \code{feature.ids}, \code{feature.labels}
#'   \item \code{cells}: one row per cell with computability metadata
#'   \item \code{landmarks}: one row per computed landmark point
#' }
#'
#' @export
linf.landmarks <- function(M,
                           csts,
                           depth = NULL,
                           view = c("active", "rare", "absorb"),
                           landmark.types = c("endpoint.max",
                                              "endpoint.min",
                                              "mean.rep",
                                              "median.rep"),
                           tie.method = c("first", "random", "error")) {
  validate.linf.csts(csts)

  view <- match.arg(view)
  tie.method <- match.arg(tie.method)
  landmark.types <- match.arg(landmark.types, several.ok = TRUE)

  X <- as.matrix(M)
  storage.mode(X) <- "numeric"

  if (any(!is.finite(X))) stop("linf.landmarks: non-finite entries found")
  if (any(X < 0, na.rm = TRUE)) stop("linf.landmarks: negative entries found")
  if (nrow(X) == 0L || ncol(X) == 0L) {
    stop("linf.landmarks: matrix has zero rows or columns")
  }
  if (nrow(X) != length(csts$cell.label)) {
    stop("linf.landmarks: nrow(M) must match the number of samples in csts")
  }

  max.depth <- csts$cst.depth %||% 1L
  if (is.null(depth)) depth <- max.depth

  if (!is.numeric(depth) || length(depth) != 1L || depth < 1L || depth %% 1 != 0) {
    stop("linf.landmarks: depth must be an integer >= 1")
  }
  if (depth > max.depth) {
    stop("linf.landmarks: requested depth exceeds csts$cst.depth")
  }

  resolved.view <- resolve.linf.landmark.view(csts, view)
  id.levels <- resolve.linf.cst.levels(csts, kind = "id", view = resolved.view)
  label.levels <- resolve.linf.cst.levels(csts, kind = "label", view = resolved.view)

  if (length(id.levels) < depth || length(label.levels) < depth) {
    stop("linf.landmarks: requested view does not contain the requested depth")
  }

  cell.ids <- as.character(id.levels[[depth]])
  cell.labels <- as.character(label.levels[[depth]])

  if (length(cell.ids) != nrow(X) || length(cell.labels) != nrow(X)) {
    stop("linf.landmarks: CST levels at the requested depth must align with rows of M")
  }

  meta <- resolve.linf.feature.meta(
    X,
    feature.ids = csts$feature.ids,
    feature.labels = csts$feature.labels
  )

  sep <- csts$sep %||% "__"
  rare.label <- csts$rare.label %||% "RARE_DOMINANT"
  row.ids <- rownames(X)
  if (is.null(row.ids)) row.ids <- rep(NA_character_, nrow(X))

  cells <- empty.linf.landmark.cells()
  landmarks <- empty.linf.landmark.rows()
  unique.cells <- unique(cell.ids[!is.na(cell.ids)])

  for (cell in unique.cells) {
    members <- which(cell.ids == cell)
    cell.label <- cell.labels[members[1L]]
    leaf.id <- tail(strsplit(cell, sep, fixed = TRUE)[[1L]], 1L)
    target <- resolve.linf.landmark.feature(
      leaf.id,
      feature.ids = meta$feature.ids,
      feature.labels = meta$feature.labels,
      rare.label = rare.label
    )

    cells <- rbind(cells, data.frame(
      cell.id = cell,
      cell.label = cell.label,
      cell.size = length(members),
      target.feature.id = target$feature.id,
      target.feature.label = target$feature.label,
      is.rare = target$is.rare,
      landmarks.computable = target$computable,
      stringsAsFactors = FALSE
    ))

    if (!isTRUE(target$computable)) next

    vals <- X[members, target$index]

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
        cell.id = cell,
        cell.label = cell.label,
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
      cells = cells,
      landmarks = landmarks
    ),
    class = "linf.landmarks"
  )
}
