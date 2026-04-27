#!/usr/bin/env Rscript

script_args <- commandArgs(trailingOnly = FALSE)
script_matches <- grep("^--file=", script_args, value = TRUE)
script_arg <- if (length(script_matches)) sub("^--file=", "", script_matches[1]) else ""
script_dir <- if (nzchar(script_arg)) {
  dirname(normalizePath(script_arg, winslash = "/", mustWork = TRUE))
} else {
  normalizePath("papers/gut-dcst/scripts", winslash = "/", mustWork = TRUE)
}

source(file.path(script_dir, "build_ibd_reviewer_sensitivity_assets.R"))
