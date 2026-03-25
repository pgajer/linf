## ============================================================================
## build_valencia2k.R
## ============================================================================
##
## Creates the bundled `valencia2k` example dataset for the linf package
## from the canonical Valencia 13k training set.
##
## Input files (not shipped with the package):
##   ~/current_projects/valencia/tx.13k.rds                          -- 12881 x 178 compositional matrix
##   ~/current_projects/valencia/cst.tx.13k.rds                      -- 12881 x 12 CST annotation data frame
##   ~/current_projects/valencia/tx.13k_sample_filtering_summary.csv -- per-sample read counts
##
## Output (saved to data/):
##   valencia2k   -- list with components:
##     $rel       -- 2000 x 178 numeric matrix (compositional, rows sum to 1)
##     $cst       -- data frame with sample_id, Val_CST, Val_subCST
##     $reads     -- integer vector of reads_after_taxa_filter (for count reconstruction)
##     $taxa      -- character vector of 178 taxon names (column names of $rel)
##     $source    -- provenance string
##
## Run from the package root:
##   source("data-raw/build_valencia2k.R")
## ============================================================================

## --- Paths ---------------------------------------------------------------

valencia_root <- path.expand("~/current_projects/valencia")

tx_path    <- file.path(valencia_root, "tx.13k.rds")
cst_path   <- file.path(valencia_root, "cst.tx.13k.rds")
sf_path    <- file.path(valencia_root, "tx.13k_sample_filtering_summary.csv")

stopifnot(file.exists(tx_path), file.exists(cst_path), file.exists(sf_path))

## --- Load ----------------------------------------------------------------

tx      <- readRDS(tx_path)                          # 12881 x 178 compositional matrix
cst_df  <- readRDS(cst_path)                         # 12881 x 12 data frame
sf      <- read.csv(sf_path, stringsAsFactors = FALSE)

## --- Validate alignment --------------------------------------------------

cst_df$sample_id <- as.character(cst_df$sample_id)
sf$sample_id     <- as.character(sf$sample_id)

stopifnot(identical(cst_df$sample_id, rownames(tx)))
stopifnot(nrow(tx) == 12881L, ncol(tx) == 178L)

## --- Read counts for the final 12881 samples ----------------------------

sf_kept <- sf[sf$kept_in_tx13k == TRUE, ]
sf_kept <- sf_kept[match(cst_df$sample_id, sf_kept$sample_id), ]
stopifnot(identical(sf_kept$sample_id, cst_df$sample_id))

## --- Stratified subsample ------------------------------------------------
##
## Draw 2000 samples stratified by Val_subCST to preserve proportional
## representation of all 13 Valencia sub-CSTs.

set.seed(42L)
n_target <- 2000L

subcst_levels <- unique(cst_df$Val_subCST)
sampled_rows  <- integer(0)

for (sub in subcst_levels) {
  idx   <- which(cst_df$Val_subCST == sub)
  frac  <- length(idx) / nrow(cst_df)
  n_sub <- max(5L, round(frac * n_target))
  n_sub <- min(n_sub, length(idx))
  sampled_rows <- c(sampled_rows, sample(idx, n_sub))
}

## Sort by original row order for reproducibility
sampled_rows <- sort(sampled_rows)

cat("Subsampled", length(sampled_rows), "samples from", nrow(tx), "\n")

## --- Build the bundle ----------------------------------------------------

rel   <- tx[sampled_rows, ]
cst_sub <- cst_df[sampled_rows, c("sample_id", "Val_CST", "Val_subCST")]
reads <- as.integer(sf_kept$reads_after_taxa_filter[sampled_rows])

## Sanity checks
stopifnot(nrow(rel) == length(sampled_rows))
stopifnot(nrow(cst_sub) == length(sampled_rows))
stopifnot(length(reads) == length(sampled_rows))
stopifnot(all(abs(rowSums(rel) - 1) < 1e-8))
stopifnot(identical(rownames(rel), cst_sub$sample_id))

## Reset row names to simple integers for cleanliness
rownames(rel)    <- paste0("s", seq_len(nrow(rel)))
cst_sub$sample_id <- rownames(rel)
rownames(cst_sub) <- NULL

## Print sub-CST distribution
cat("\nVal_subCST distribution in subsample:\n")
print(sort(table(cst_sub$Val_subCST), decreasing = TRUE))

cat("\nVal_CST distribution in subsample:\n")
print(sort(table(cst_sub$Val_CST), decreasing = TRUE))

## --- Package as a single list -------------------------------------------

valencia2k <- list(
  rel    = rel,
  cst    = cst_sub,
  reads  = reads,
  taxa   = colnames(rel),
  source = paste0(
    "Stratified subsample (n=", nrow(rel),
    ") from the Valencia 13k vaginal 16S training set ",
    "(France et al. 2020, doi:10.1128/mSystems.00149-20). ",
    "Original dataset: 12881 samples x 178 taxa. ",
    "Subsampled with set.seed(42), stratified by Val_subCST. ",
    "Matrix is compositional (rows sum to 1). ",
    "Use sweep(valencia2k$rel, 1, valencia2k$reads, '*') to reconstruct count-like values."
  )
)

## --- Save ----------------------------------------------------------------

save(valencia2k, file = "data/valencia2k.rda", compress = "xz")

cat("\nSaved data/valencia2k.rda\n")
cat("File size:", round(file.info("data/valencia2k.rda")$size / 1024), "KB\n")
