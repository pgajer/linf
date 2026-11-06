source("tools/qa-paths.R")
output <- qa_dir("site")
pkgdown::build_site(preview = FALSE, new_process = TRUE,
                   override = list(destination = output))
message("Website preview: ", file.path(output, "index.html"))
