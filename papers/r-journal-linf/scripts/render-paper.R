#!/usr/bin/env Rscript

required <- c("rmarkdown", "rjtools", "linf", "Matrix", "knitr", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  stop("Missing build package(s): ", paste(missing, collapse = ", "))
}

target <- package_version(jsonlite::fromJSON("package-source.json")$version)
if (utils::packageVersion("linf") != target) {
  stop("The article requires linf ", target, "; found ", utils::packageVersion("linf"))
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

# The combined format's internal PDF pass cleans its supporting figures.
# Render explicitly with cleanup disabled so the TeX sources remain usable.
pdf_environment <- new.env(parent = globalenv())
rmarkdown::render(
  input = "linf.Rmd",
  output_format = "rjtools::rjournal_pdf_article",
  envir = pdf_environment,
  clean = FALSE,
  quiet = FALSE
)

runtime <- utils::sessionInfo()
machine <- as.list(Sys.info()[c("sysname", "release", "machine")])
machine$logical_cores <- parallel::detectCores(logical = TRUE)
if (identical(Sys.info()[["sysname"]], "Darwin")) {
  machine$cpu_model <- unname(system2("sysctl", c("-n", "machdep.cpu.brand_string"),
                                    stdout = TRUE))
  machine$memory_bytes <- as.numeric(system2("sysctl", c("-n", "hw.memsize"),
                                           stdout = TRUE))
}
versions <- vapply(sort(loadedNamespaces()), function(package) {
  utils::packageDescription(package, fields = "Version")
}, character(1L))
utils::write.csv(pdf_environment$benchmark_results,
                 "build/benchmark-results.csv", row.names = FALSE)
writeLines(c("Computational session for the PDF benchmark results", "",
             capture.output(print(runtime))), "build/session-info.txt")
jsonlite::write_json(list(
  recorded_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  artifact = "linf.pdf",
  results = "build/benchmark-results.csv",
  R = R.version.string,
  platform = R.version$platform,
  machine = machine,
  packages = as.list(versions)
), "build/benchmark-environment.json", pretty = TRUE, auto_unbox = TRUE)

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
    paste("Matrix version:", utils::packageDescription("Matrix", fields = "Version")),
    sprintf("render seconds: %.3f", elapsed)
  ),
  "build/render-info.txt"
)

message(sprintf("Article render completed in %.2f seconds.", elapsed))
