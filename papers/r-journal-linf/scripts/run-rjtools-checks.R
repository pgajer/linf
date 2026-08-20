#!/usr/bin/env Rscript

if (!requireNamespace("rjtools", quietly = TRUE)) {
  stop("rjtools is required")
}

dir.create("build", showWarnings = FALSE, recursive = TRUE)
logfile <- normalizePath("build", mustWork = TRUE)
logfile <- file.path(logfile, "rjtools-initial-checks.log")
if (file.exists(logfile)) file.remove(logfile)

stage <- tempfile("linf-rjtools-check-")
dir.create(stage, recursive = TRUE)
on.exit(unlink(stage, recursive = TRUE, force = TRUE), add = TRUE)

top_level <- c(
  "linf.Rmd", "linf.tex", "linf.R", "linf.pdf", "RJreferences.bib", "RJournal.sty",
  "_Rpackages.txt"
)
if (!all(file.copy(top_level, stage, overwrite = TRUE))) {
  stop("Could not stage the generated article files for rjtools checks")
}
dir.create(file.path(stage, "motivation-letter"))
letter_files <- c(
  "motivation-letter/motivation-letter.md",
  "motivation-letter/motivation-letter.pdf"
)
if (!all(file.copy(
  letter_files,
  file.path(stage, "motivation-letter"),
  overwrite = TRUE
))) {
  stop("Could not stage the motivation letter for rjtools checks")
}

# rjtools 1.0.21 passes the base `file` function to check_date() from
# initial_check_article(). Run the same checks explicitly so that the date
# checker receives the generated article file it expects.
log_connection <- file(logfile, open = "a", encoding = "UTF-8")
on.exit(close(log_connection), add = TRUE)

old_options <- options(
  check.log.file = log_connection,
  check.log.journal = new.env(parent = emptyenv())
)
on.exit(options(old_options), add = TRUE)

checks <- list(
  filenames = rjtools::check_filenames(stage),
  structure = rjtools::check_structure(stage),
  folders = rjtools::check_folder_structure(stage),
  unnecessary_files = rjtools::check_unnecessary_files(stage),
  cover_letter = rjtools::check_cover_letter(stage),
  title = rjtools::check_title(stage),
  sections = rjtools::check_section(stage),
  abstract = rjtools::check_abstract(stage),
  spelling = rjtools::check_spelling(
    stage,
    ignore = c(
      "linf", "compositional", "microbiome", "dcst", "asv", "csts",
      "dcsts", "iter", "valencia", "agp", "rrrrrr", "densify", "asize",
      "nonnegative", "interpretable", "centroid", "centroids", "backend",
      "backends", "thresholded", "dataset", "datasets"
    )
  ),
  proposed_package = rjtools::check_proposed_pkg("linf", ask = FALSE),
  package_labels = rjtools::check_pkg_label(stage),
  package_availability = rjtools::check_packages_available(stage),
  bibliography = rjtools::check_bib_doi(stage),
  csl = rjtools::check_csl(stage),
  date = rjtools::check_date(stage, file.path(stage, "linf.tex"))
)

result <- getOption("check.log.journal")$results
capture.output(str(list(return_values = checks, results = result)),
               file = "build/rjtools-results.txt")

flattened <- unlist(result, recursive = TRUE, use.names = TRUE)
error_values <- grepl("ERROR|FAIL", as.character(flattened), ignore.case = TRUE)
if (any(error_values)) {
  stop(
    "rjtools reported an error; see ", logfile, " and build/rjtools-results.txt"
  )
}

message("All applicable rjtools checks, including the corrected date check, completed without an error status.")
