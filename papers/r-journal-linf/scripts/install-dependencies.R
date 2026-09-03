#!/usr/bin/env Rscript

packages <- setdiff(readLines("_Rpackages.txt", warn = FALSE), "linf")
packages <- packages[nzchar(packages)]
library <- Sys.getenv("R_LIBS_USER")
if (!nzchar(library)) stop("Set R_LIBS_USER to the article library")
dir.create(library, recursive = TRUE, showWarnings = FALSE)
.libPaths(c(library, .libPaths()))
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  install.packages(missing, lib = library, repos = "https://cloud.r-project.org")
}
stopifnot(all(vapply(packages, requireNamespace, logical(1), quietly = TRUE)))
message("Article dependencies are available.")
