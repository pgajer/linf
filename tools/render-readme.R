source("tools/qa-paths.R")
pkgload::load_all(".", quiet = TRUE)
# README.md is maintained output; temporary files belong outside the source tree.
rmarkdown::render("README.Rmd", output_options = list(html_preview = FALSE),
                  intermediates_dir = qa_dir("readme"), quiet = TRUE)
