# Independent graph and rendered-trace checks for the website-only article.
# Run from the package root; package checks exclude this optional article.
source <- "vignettes/articles/valencia-hypercube-embedding.Rmd"
script <- tempfile(fileext = ".R")
invisible(knitr::purl(source, output = script, documentation = 0L, quiet = TRUE))
expressions <- parse(script)
unlink(script)
e <- new.env()
for (expr in expressions) {
  if (is.call(expr) && identical(expr[[1]], as.name("<-")) && is.symbol(expr[[2]]) &&
      as.character(expr[[2]]) %in% c("knn_edges", "graph_components", "plot_embedding")) {
    eval(expr, e)
  }
}
stopifnot(all(c("knn_edges", "graph_components", "plot_embedding") %in% ls(e)))
canonical <- function(edges) {
  x <- t(apply(edges, 1L, sort))
  unname(x[order(x[, 1], x[, 2]), , drop = FALSE])
}
X <- matrix(c(0, 10, 11.1, 12.4), ncol = 1)
g <- e$knn_edges(X, 1)
stopifnot(identical(g$edges, cbind(1:3, 2:4)), length(unique(e$graph_components(4, g$edges))) == 1L)
p <- 4:1
stopifnot(identical(canonical(matrix(p[e$knn_edges(X[p, , drop = FALSE], 1)$edges], ncol = 2)),
                    canonical(g$edges)))
# Equal-distance neighbors must all survive reordering, and duplicates have
# positive optimizer weights without creating self edges.
tied <- matrix(c(0, 1, -1, 0), ncol = 1)
a <- e$knn_edges(tied, 1)
b <- e$knn_edges(tied[p, , drop = FALSE], 1)
stopifnot(identical(canonical(a$edges), canonical(matrix(p[b$edges], ncol = 2))),
          all(a$edges[, 1] < a$edges[, 2]), all(a$weights > 0),
          nrow(e$knn_edges(matrix(1, 1, 1), 1)$edges) == 0L,
          length(unique(e$graph_components(4, matrix(c(1L, 3L, 2L, 4L), ncol = 2)))) == 2L)
pkgload::load_all(".", quiet = TRUE)
e$meta <- valencia_linf_hypercube_1k$meta
e$component.colors <- c(Li = "#0072B2", Lc = "#E69F00", Gv = "#009E73", Bv = "#CC79A7")
X4 <- valencia_linf_hypercube_1k$rel4
for (ref in colnames(X4)) {
  coords <- linf.hypercube.embedding(X4, ref)
  widget <- e$plot_embedding(coords, ref, "audit", unit.cube = TRUE)
  built <- withCallingHandlers(plotly::plotly_build(widget), warning = function(w) stop(w))
  traces <- built$x$data
  stopifnot(sum(vapply(traces, function(x) length(x$x), integer(1))) == nrow(X4),
            length(traces) == 4L, all(vapply(traces, function(x) identical(x$marker$symbol, "circle"), logical(1))))
  hover <- unlist(lapply(traces, `[[`, "text"))
  plotted.ids <- sub("<br>.*$", "", sub("^sample: ", "", hover))
  stopifnot(identical(sort(unname(plotted.ids)), sort(as.character(e$meta$sample_id))))
  for (category in unique(e$meta$Val_CST)) {
    stopifnot(any(grepl(paste0("VALENCIA CST: ", category, "<br>"), hover, fixed = TRUE)))
  }
}
cat("Union graph: asymmetric relations, ties, reordering, duplicates and components verified.\n")
cat("Four native widgets: all 1000 points, four component groups, seven CST categories; no build warnings.\n")
