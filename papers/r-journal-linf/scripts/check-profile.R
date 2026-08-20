tidy_dir <- Sys.getenv("LINF_TIDY_DIR", unset = "")
if (nzchar(tidy_dir) && file.exists(file.path(tidy_dir, "tidy"))) {
  Sys.setenv(PATH = paste(tidy_dir, Sys.getenv("PATH"), sep = .Platform$path.sep))
}
