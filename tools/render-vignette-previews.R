# Render installed-vignette sources into an excluded local preview directory.
# No source data downloads or interactive article dependencies are needed.
pkgload::load_all(".", quiet = TRUE)
output <- normalizePath("build", mustWork = FALSE)
dir.create(file.path(output, "vignettes"), recursive = TRUE, showWarnings = FALSE)
for (source in list.files("vignettes", "[.]Rmd$", full.names = TRUE)) {
  rmarkdown::render(source, output_dir = file.path(output, "vignettes"),
                    envir = new.env(parent = globalenv()), quiet = TRUE)
}
