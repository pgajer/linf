# Render installed-vignette sources outside the package source tree.
# No source data downloads or interactive article dependencies are needed.
pkgload::load_all(".", quiet = TRUE)
source("tools/qa-paths.R")
output <- qa_dir("vignettes")
for (source in list.files("vignettes", "[.]Rmd$", full.names = TRUE)) {
  rmarkdown::render(source, output_dir = output,
                    envir = new.env(parent = globalenv()), quiet = TRUE)
}
