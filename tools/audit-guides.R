# Run from the package root: make audit-guides
# Check maintained catalogs, source definitions and installed help coverage.
ns <- readLines("NAMESPACE", warn = FALSE)
exports <- sub("^export[(](.*)[)]$", "\\1", grep("^export[(]", ns, value = TRUE))
methods <- sub("^S3method[(]([^,]+),([^)]*)[)]$", "\\1.\\2",
               grep("^S3method[(]", ns, value = TRUE))
guide <- readLines("vignettes/function-guide.Rmd", warn = FALSE)
catalog <- sub("^\\| `([^`]+)[(][)]`.*$", "\\1",
               grep("^\\| `[^`]+[(][)]` \\|", guide, value = TRUE))
stopifnot(!anyDuplicated(catalog), setequal(catalog, exports))
stopifnot(any(grepl(sprintf("%d explicit function exports", length(exports)),
                   guide, fixed = TRUE)),
          length(methods) == 2L, length(intersect(exports, methods)) == 0L)
definitions <- character()
for (path in list.files("R", pattern = "[.]R$", full.names = TRUE)) {
  for (expr in parse(path)) {
    if (is.call(expr) && identical(expr[[1L]], as.name("<-")) &&
        length(expr) == 3L && is.call(expr[[3L]]) &&
        identical(expr[[3L]][[1L]], as.name("function"))) {
      definitions <- c(definitions, as.character(expr[[2L]]))
    }
  }
}
stopifnot(all(c(exports, methods) %in% definitions))
# Methods that are also exports still need exactly one ordinary catalog row.
stopifnot(all(intersect(exports, methods) %in% catalog))
registered.generics <- sub("[.]linf[.]csts$", "", methods)
stopifnot(all(vapply(registered.generics, function(generic) {
  any(grepl(paste0("`", generic, "()"), guide, fixed = TRUE))
}, logical(1))))

rd <- lapply(list.files("man", "[.]Rd$", full.names = TRUE), tools::parse_Rd)
aliases <- unlist(lapply(rd, function(doc) {
  unlist(lapply(doc, function(node) {
    if (identical(attr(node, "Rd_tag"), "\\alias")) paste(unlist(node), collapse = "")
  }))
}))
data.names <- unlist(lapply(list.files("data", "[.]rda$", full.names = TRUE),
                           function(path) load(path, envir = new.env())))
datasets <- readLines("vignettes/example-datasets.Rmd", warn = FALSE)
data.catalog <- sub("^\\| `([^`]+)`.*$", "\\1",
                    grep("^\\| `[^`]+` \\|", datasets, value = TRUE))
stopifnot(!anyDuplicated(data.catalog), setequal(data.catalog, data.names),
          length(data.names) == 5L,
          all(c(exports, methods, data.names, "linf", "linf-package") %in% aliases))
for (source in c("vignettes/function-guide.Rmd", "vignettes/example-datasets.Rmd")) {
  text <- readLines(source, warn = FALSE)
  stopifnot(any(grepl("%\\VignetteIndexEntry{", text, fixed = TRUE)),
            any(grepl("%\\VignetteEngine{knitr::rmarkdown}", text, fixed = TRUE)))
  links <- regmatches(text, gregexpr("\\]\\([^():]+[.]html(#[^)]*)?\\)", text))
  links <- sub("#.*$", "", sub("^\\]\\(|\\)$", "", unlist(links)))
  links <- sub("\\)$", "", links)
  if (length(links)) stopifnot(all(file.exists(file.path("vignettes",
                                          sub("[.]html$", ".Rmd", links)))))
}
cat(sprintf("Catalog verified: %d explicit functions, %d registered methods (%d overlapping exports), %d datasets.\n",
            length(exports), length(methods), length(intersect(exports, methods)), length(data.names)))
cat("Source definitions, help aliases, vignette metadata and local guide links verified.\n")

# Alias presence alone does not make the overview discoverable in installed help.
package.rd <- readLines("man/linf-package.Rd", warn = FALSE)
stopifnot(!any(grepl("\\keyword{internal}", package.rd, fixed = TRUE)))
args <- commandArgs(trailingOnly = TRUE)
if (length(args)) {
  stopifnot(length(args) == 1L, startsWith(args, "--library="))
  lib <- sub("^--library=", "", args)
  package <- find.package("linf", lib.loc = lib)
  stopifnot(length(utils::help("linf", package = "linf", lib.loc = lib)) == 1L,
            length(utils::help("linf-package", package = "linf", lib.loc = lib)) == 1L)
  provenance <- file.path(package, "DATA_PROVENANCE.md")
  manifest <- file.path(package, "DATA_MANIFEST.csv")
  stopifnot(file.exists(provenance), file.exists(manifest),
            any(grepl("Reproduction record", readLines(provenance), fixed = TRUE)),
            nrow(read.csv(manifest)) >= 10L)
  index <- readLines(file.path(package, "html", "00Index.html"), warn = FALSE)
  stopifnot(any(grepl("linf-package.html", index, fixed = TRUE)))
  items <- utils::vignette(package = "linf", lib.loc = lib)$results[, "Item"]
  expected <- c("function-guide", "example-datasets", "linf-intro", "linf-vaginal")
  stopifnot(setequal(items, expected))
  for (item in expected) {
    html <- file.path(package, "doc", paste0(item, ".html"))
    stopifnot(file.exists(html))
    text <- readLines(html, warn = FALSE)
    links <- regmatches(text, gregexpr('href="[^":]+[.]html(#[^"]*)?"', text))
    links <- sub('^href="|"$', "", unlist(links))
    links <- sub('"$', "", sub("#.*$", "", links))
    local <- links[basename(links) %in% paste0(expected, ".html")]
    stopifnot(all(file.exists(file.path(package, "doc", local))))
  }
  cat("Installed package overview, both aliases, four vignettes and local guide links verified.\n")
}
