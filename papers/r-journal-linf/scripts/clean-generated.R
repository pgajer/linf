#!/usr/bin/env Rscript

files <- c(
  "linf.html", "linf.pdf", "linf.tex", "linf.R",
  "RJournal.tex", "RJournal.pdf", "Rplots.pdf",
  "RJwrapper.tex", "RJwrapper.pdf",
  "linf.log", "linf.aux", "linf.out", "linf.bbl", "linf.blg",
  "RJwrapper.log", "RJwrapper.aux", "RJwrapper.out",
  "motivation-letter/motivation-letter.pdf"
)
dirs <- c("linf_files", "build", "output", "library")
unlink(files[file.exists(files)], force = TRUE)
unlink(dirs[dir.exists(dirs)], recursive = TRUE, force = TRUE)
message("Removed generated article artifacts.")
