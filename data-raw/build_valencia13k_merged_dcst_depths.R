## ============================================================================
## build_valencia13k_merged_dcst_depths.R
## ============================================================================
##
## Creates lightweight bundled assets for depth-2 and depth-3 merged dCST
## assignments on the canonical Valencia 13k training set.
##
## Input files (not shipped with the package):
##   $LINF_VALENCIA_ROOT/tx.13k.rds      -- compositional matrix
##   $LINF_VALENCIA_ROOT/cst.tx.13k.rds  -- CST annotation data frame
##
## Output (saved to data/):
##   valencia13k_dcst_depth2_merged
##   valencia13k_dcst_depth3_merged
##
## Run from the package root:
##   source("data-raw/build_valencia13k_merged_dcst_depths.R")
## ============================================================================

if (!exists("linf.csts", mode = "function")) {
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("pkgload is required when running this script from the source tree")
  }
  pkgload::load_all(".", quiet = TRUE)
}

valencia_root <- Sys.getenv("LINF_VALENCIA_ROOT")
if (!nzchar(valencia_root)) {
  stop("Set LINF_VALENCIA_ROOT to the directory containing the Valencia source files.")
}
valencia_root <- normalizePath(valencia_root, mustWork = TRUE)

tx_path <- file.path(valencia_root, "tx.13k.rds")
cst_path <- file.path(valencia_root, "cst.tx.13k.rds")

stopifnot(file.exists(tx_path), file.exists(cst_path))

tx <- readRDS(tx_path)
cst_df <- readRDS(cst_path)

cst_df$sample_id <- as.character(cst_df$sample_id)
stopifnot(identical(cst_df$sample_id, rownames(tx)))

source_n <- nrow(tx)
source_p <- ncol(tx)
anon_sample_id <- sprintf("v13k_%05d", seq_len(source_n))

rownames(tx) <- anon_sample_id

n0 <- 50L
refinement.factor <- 2
sep <- "__"
rare.label <- "RARE_DOMINANT"

linf.rel <- normalize.linf(tx)

dcst1 <- linf.csts(
  linf.rel,
  n0 = n0,
  low.freq.policy = "absorb",
  rare.label = rare.label
)

dcst2 <- refine.linf.csts(
  linf.rel,
  dcst1,
  n0 = n0,
  refinement.factor = refinement.factor,
  sep = sep,
  low.freq.policy = "absorb",
  rare.label = rare.label,
  verbose = FALSE
)

dcst3 <- refine.linf.csts.iter(
  linf.rel,
  dcst2,
  n0 = n0,
  refinement.factor = refinement.factor,
  sep = sep,
  low.freq.policy = "absorb",
  rare.label = rare.label,
  verbose = FALSE
)

summarize_labels <- function(labels, depth, sep) {
  tab <- sort(table(labels), decreasing = TRUE)
  data.frame(
    depth = as.integer(depth),
    dcst_label = names(tab),
    n = as.integer(tab),
    prop = as.numeric(tab) / sum(tab),
    path_length = vapply(
      strsplit(names(tab), sep, fixed = TRUE),
      length,
      integer(1)
    ),
    stringsAsFactors = FALSE,
    row.names = NULL
  )
}

make_assignments <- function(csts) {
  depth <- csts$depth
  out <- data.frame(
    sample_id = anon_sample_id,
    source_row = seq_len(source_n),
    Val_CST = as.character(cst_df$Val_CST),
    Val_subCST = as.character(cst_df$Val_subCST),
    stringsAsFactors = FALSE
  )

  for (d in seq_len(depth)) {
    out[[paste0("dcst_depth", d)]] <- as.character(csts$lineage.labels.absorb[[d]])
  }

  out
}

make_asset <- function(csts, depth) {
  assignments <- make_assignments(csts)
  summaries <- setNames(
    lapply(seq_len(depth), function(d) {
      summarize_labels(csts$lineage.labels.absorb[[d]], depth = d, sep = sep)
    }),
    paste0("depth", seq_len(depth))
  )

  list(
    assignments = assignments,
    summaries = summaries,
    feature_labels = colnames(tx),
    params = list(
      source_n = source_n,
      source_p = source_p,
      n0 = n0,
      low.freq.policy = "absorb",
      refinement.factor = refinement.factor,
      sep = sep,
      rare.label = rare.label,
      sample_id_policy = "anonymized sequential IDs; source_row records source row number"
    ),
    source = paste0(
      "Depth-", depth,
      " merged dCST assignment asset generated from the Valencia 13k vaginal ",
      "16S training set. Input matrix: tx.13k.rds with ", source_n,
      " samples x ", source_p, " taxa. L-infinity normalized with ",
      "normalize.linf(), then computed with n0=50, low.freq.policy='absorb', ",
      "and refinement.factor=2 for hierarchical refinements. Original source ",
      "sample IDs are not stored in the bundled object."
    )
  )
}

valencia13k_dcst_depth2_merged <- make_asset(dcst2, depth = 2L)
valencia13k_dcst_depth3_merged <- make_asset(dcst3, depth = 3L)

stopifnot(nrow(valencia13k_dcst_depth2_merged$assignments) == source_n)
stopifnot(nrow(valencia13k_dcst_depth3_merged$assignments) == source_n)
stopifnot(!any(grepl(rare.label, valencia13k_dcst_depth2_merged$assignments$dcst_depth2, fixed = TRUE)))
stopifnot(!any(grepl(rare.label, valencia13k_dcst_depth3_merged$assignments$dcst_depth3, fixed = TRUE)))

cat("Depth-2 merged dCST count:", nrow(valencia13k_dcst_depth2_merged$summaries$depth2), "\n")
cat("Depth-3 merged dCST count:", nrow(valencia13k_dcst_depth3_merged$summaries$depth3), "\n")

save(
  valencia13k_dcst_depth2_merged,
  file = "data/valencia13k_dcst_depth2_merged.rda",
  compress = "xz"
)

save(
  valencia13k_dcst_depth3_merged,
  file = "data/valencia13k_dcst_depth3_merged.rda",
  compress = "xz"
)

cat("\nSaved data/valencia13k_dcst_depth2_merged.rda\n")
cat("File size:", round(file.info("data/valencia13k_dcst_depth2_merged.rda")$size / 1024), "KB\n")
cat("Saved data/valencia13k_dcst_depth3_merged.rda\n")
cat("File size:", round(file.info("data/valencia13k_dcst_depth3_merged.rda")$size / 1024), "KB\n")
