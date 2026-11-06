# Shared output routing for developer tools; never used by installed package code.
qa_dir <- function(..., create = TRUE) {
  root <- Sys.getenv("LINF_QA_DIR", unset = "")
  if (!nzchar(root) && nzchar(Sys.which("git"))) {
    configured <- suppressWarnings(system2(
      "git", c("config", "--get", "linf.qaDir"), stdout = TRUE, stderr = FALSE
    ))
    if (length(configured) == 1L) root <- configured
  }
  if (!nzchar(root)) root <- file.path(tools::R_user_dir("linf", "cache"), "qa")
  root <- path.expand(root)
  if (!grepl("^(/|[A-Za-z]:[/\\\\])", root)) root <- file.path(getwd(), root)
  root <- normalizePath(root, winslash = "/", mustWork = FALSE)
  source_root <- normalizePath(getwd(), winslash = "/", mustWork = TRUE)
  if (identical(root, source_root) || startsWith(root, paste0(source_root, "/"))) {
    stop("Set the QA output directory outside the package source tree.")
  }
  path <- file.path(root, ...)
  if (create) dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

if (identical(commandArgs(trailingOnly = TRUE), "--print")) {
  cat(qa_dir(create = FALSE), "\n", sep = "")
}
