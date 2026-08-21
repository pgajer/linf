#!/usr/bin/env Rscript

required <- c("rmarkdown", "rjtools", "linf", "Matrix", "knitr")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing build package(s): ", paste(missing, collapse = ", "))
}

if (utils::packageVersion("linf") < "0.3.0") {
  stop("The article requires linf >= 0.3.0; found ", utils::packageVersion("linf"))
}

dir.create("build", showWarnings = FALSE, recursive = TRUE)
started <- proc.time()[["elapsed"]]

rmarkdown::render(
  input = "linf.Rmd",
  output_format = "rjtools::rjournal_article",
  envir = new.env(parent = globalenv()),
  clean = FALSE,
  quiet = FALSE
)

rmarkdown::render(
  input = "motivation-letter/motivation-letter.md",
  output_file = "motivation-letter.pdf",
  output_dir = "motivation-letter",
  clean = TRUE,
  quiet = TRUE
)

required_outputs <- c(
  "linf.html",
  "linf.pdf",
  "linf.tex",
  "linf.R",
  "motivation-letter/motivation-letter.pdf"
)
missing_outputs <- required_outputs[!file.exists(required_outputs)]
if (length(missing_outputs)) {
  stop("Article build did not create: ", paste(missing_outputs, collapse = ", "))
}

elapsed <- proc.time()[["elapsed"]] - started
writeLines(
  c(
    paste("linf version:", as.character(utils::packageVersion("linf"))),
    paste("R version:", R.version.string),
    paste("rjtools version:", as.character(utils::packageVersion("rjtools"))),
    sprintf("render seconds: %.3f", elapsed)
  ),
  "build/render-info.txt"
)

message(sprintf("Article render completed in %.2f seconds.", elapsed))
