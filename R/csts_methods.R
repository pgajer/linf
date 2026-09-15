#' Print method for \code{"linf.csts"}
#'
#' @param x A \code{"linf.csts"} object.
#' @param ... Unused.
#' @return The input object, invisibly. Printing shows sample accounting and up
#'   to ten groups per depth. Use \code{table(x$lineage.labels[[d]],
#'   useNA = "ifany")} for the complete counts at depth \code{d}.
#' @examples
#' M <- rbind(
#'   s1 = c(A = 1.0, B = 0.2),
#'   s2 = c(A = 0.8, B = 1.0)
#' )
#' fit <- linf.csts(M, n0 = 1)
#' print(fit)
#' @export
print.linf.csts <- function(x, ...) {

  validate.linf.csts(x)

  original <- x
  x <- linf.ensure.nodes(x)
  stats <- summary(x)
  cat("Dominant Community State Type Hierarchy\n")
  cat("Total samples:", length(x$lineage.label), " | Depth:", x$depth,
      " | Policy:", x$low.freq.policy, "\n")
  for (d in seq_len(x$depth)) {
    cat("\nDepth", d, "-", stats$assigned.samples[d], "assigned,",
        stats$unassigned.samples[d], "unassigned,", stats$rare.samples[d], "rare\n")
    tab <- sort(table(x$lineage.labels[[d]]), decreasing = TRUE)
    if (!length(tab)) {
      cat("  No assigned lineages at this depth.\n")
    } else {
      if (stats$retained.lineages[d] == 0L) cat("  No retained lineages; all assignments are synthetic rare groups.\n")
      for (nm in utils::head(names(tab), 10L)) cat("  ", nm, ": ", tab[[nm]], "\n", sep = "")
      if (length(tab) > 10L) {
        cat(sprintf("  ... %d more lineages. Use table(x$lineage.labels[[%d]], useNA = 'ifany') for all counts.\n",
                    length(tab) - 10L, d))
      }
    }
  }

  invisible(original)
}

#' Summarize dCST hierarchy statistics
#'
#' @param object A \code{"linf.csts"} object.
#' @param ... Unused.
#' @return Data frame with one row per stored depth. \code{input.samples},
#'   \code{assigned.samples}, \code{unassigned.samples} and \code{rare.samples}
#'   expose the denominators. \code{n.lineages} includes synthetic rare groups;
#'   \code{retained.lineages} excludes them. For compatibility,
#'   \code{total.samples} remains an alias of \code{assigned.samples}.
#'   Size statistics include rare groups and are \code{NA} when no assignments
#'   exist; valid empty fits do not produce warnings.
#' @examples
#' M <- rbind(
#'   s1 = c(A = 1.0, B = 0.2),
#'   s2 = c(A = 0.9, B = 0.3),
#'   s3 = c(A = 0.2, B = 1.0)
#' )
#' fit <- linf.csts(M, n0 = 1)
#' summary(fit)
#' @export
summary.linf.csts <- function(object, ...) {

  validate.linf.csts(object)

  object <- linf.ensure.nodes(object)
  rows <- lapply(seq_len(object$depth), function(d) {
    keys <- object$node.ids[[d]]
    nodes <- object$nodes[[d]]
    tab <- table(keys)
    assigned <- sum(!is.na(keys))
    rare <- sum(keys %in% nodes$node.id[nodes$is.rare])
    data.frame(
      depth = as.integer(d),
      n.lineages = length(tab),
      total.samples = assigned,
      mean.size = if (length(tab)) mean(tab) else NA_real_,
      median.size = if (length(tab)) stats::median(tab) else NA_real_,
      min.size = if (length(tab)) min(tab) else NA_real_,
      max.size = if (length(tab)) max(tab) else NA_real_,
      input.samples = length(keys),
      assigned.samples = assigned,
      unassigned.samples = sum(is.na(keys)),
      rare.samples = rare,
      retained.lineages = sum(!nodes$is.rare)
    )
  })
  do.call(rbind, rows)
}
