#!/usr/bin/env Rscript

pin <- jsonlite::fromJSON("package-source.json")
target <- package_version(pin$version)
installed <- utils::packageVersion("linf")

available <- utils::available.packages(
  repos = c(CRAN = "https://cloud.r-project.org"),
  type = "source"
)
if (!"linf" %in% rownames(available)) {
  stop("linf is not listed in the current CRAN source index")
}
public <- package_version(available["linf", "Version"])

cat("Target source version:", as.character(target), "\n")
cat("Installed article version:", as.character(installed), "\n")
cat("Public CRAN version:", as.character(public), "\n")

if (installed != target) {
  stop("The isolated article library does not contain the target source version")
}
if (public < target) {
  stop(
    "Submission blocked: CRAN has linf ", public,
    " but the article describes linf ", target
  )
}

message("Public package-version gate passed.")
