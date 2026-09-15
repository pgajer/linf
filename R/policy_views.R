#' Select a stored dCST policy view
#'
#' @description
#' Returns a copy of a \code{"linf.csts"} object using either the stored
#' \code{"pure"} or \code{"absorb"} hierarchy. This does not recompute dCSTs.
#'
#' @param csts A \code{"linf.csts"} object produced by \code{\link{linf.csts}}
#'   and optionally refined by repeated calls to
#'   \code{\link{refine.linf.csts}}.
#' @param view Character. The stored policy view to activate: \code{"pure"} or
#'   \code{"absorb"}.
#'
#' @return A \code{"linf.csts"} object using the requested view.
#'
#' @examples
#' M <- rbind(
#'   s1 = c(A = 1.0, B = 0.2, C = 0.1),
#'   s2 = c(A = 0.9, B = 0.3, C = 0.1),
#'   s3 = c(A = 0.2, B = 1.0, C = 0.1),
#'   s4 = c(A = 0.2, B = 0.1, C = 1.0)
#' )
#' fit <- linf.csts(M, n0 = 2, low.freq.policy = "pure")
#' table(fit$lineage.label)
#' table(dcst.view(fit, view = "absorb")$lineage.label)
#'
#' @export
dcst.view <- function(csts, view = c("absorb", "pure")) {

  validate.linf.csts(csts)
  csts <- linf.ensure.nodes(csts)
  view <- match.arg(view)

  lineage.label <- csts[[paste0("lineage.label.", view)]]
  lineage.id <- csts[[paste0("lineage.id.", view)]]
  depth1.feature.index <- csts[[paste0("depth1.feature.index.", view)]]
  depth1.feature.id <- csts[[paste0("depth1.feature.id.", view)]]
  depth1.feature.label <- csts[[paste0("depth1.feature.label.", view)]]
  lineage.ids <- csts[[paste0("lineage.ids.", view)]]
  lineage.labels <- csts[[paste0("lineage.labels.", view)]]

  if (is.null(lineage.label) || is.null(lineage.ids) || is.null(lineage.labels)) {
    stop(sprintf('dcst.view: object does not contain the "%s" view', view))
  }

  csts$depth1.feature.index <- depth1.feature.index
  csts$depth1.feature.id <- depth1.feature.id
  csts$depth1.feature.label <- depth1.feature.label
  csts$lineage.id <- lineage.id
  csts$lineage.label <- lineage.label
  csts$lineage.ids <- lineage.ids
  csts$lineage.labels <- lineage.labels
  csts$low.freq.policy <- view
  csts <- linf.activate.nodes(csts)

  csts$landmarks <- NULL
  class(csts) <- "linf.csts"
  csts
}
