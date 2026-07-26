## ============================================================================
## build_agp_gut.R
## ============================================================================
##
## Creates the bundled `agp_gut` example dataset for the linf package from
## the deterministic AGP demonstration-subset files.
##
## `data-raw/create_agp_gut_subset.py` deterministically selects membership
## from the retained 5,000-sample dCST assignments, then writes the count and
## metadata inputs below. Phenotypes are joined after membership is fixed.
##
## Input files (in inst/extdata/):
##   agp_gut_counts.csv.gz  -- 766 x 314 count matrix (samples x taxa)
##   agp_gut_meta.csv       -- 766-row metadata with dCST labels and disease indicators
##   agp_gut_taxa.txt       -- 314 SILVA taxonomy strings
##
## Output (saved to data/):
##   agp_gut  -- list with components:
##     $counts  -- 766 x 314 integer matrix (raw counts, rows = samples, cols = taxa)
##     $meta    -- data frame with Run, dcst_depth1, dcst_depth2, disease indicators, BMI
##     $taxa    -- character vector of 314 SILVA taxonomy strings
##     $source  -- provenance string
##
## Run from the package root:
##   source("data-raw/build_agp_gut.R")
## ============================================================================

## --- Load CSV files --------------------------------------------------------

counts_df <- read.csv("inst/extdata/agp_gut_counts.csv.gz", row.names = 1,
                       check.names = FALSE)
meta_df   <- read.csv("inst/extdata/agp_gut_meta.csv", stringsAsFactors = FALSE)
taxa      <- readLines("inst/extdata/agp_gut_taxa.txt")

## --- Validate --------------------------------------------------------------

stopifnot(nrow(counts_df) == nrow(meta_df))
stopifnot(ncol(counts_df) == length(taxa))
stopifnot(all(rownames(counts_df) == meta_df$Run))

cat("Loaded:", nrow(counts_df), "samples x", ncol(counts_df), "taxa\n")

## --- Build the count matrix ------------------------------------------------

counts <- as.matrix(counts_df)
storage.mode(counts) <- "integer"

## Sanity: no negative values, all rows have reads
stopifnot(all(counts >= 0L, na.rm = TRUE))
stopifnot(all(rowSums(counts) > 0L))

cat("Library sizes: min =", min(rowSums(counts)),
    ", median =", median(rowSums(counts)),
    ", max =", max(rowSums(counts)), "\n")

## --- Print dCST distribution -----------------------------------------------

cat("\ndCST distribution:\n")
print(sort(table(meta_df$dcst_depth1), decreasing = TRUE))

cat("\nDisease prevalence:\n")
disease_cols <- c("IBS", "IBD", "Diabetes", "Autoimmune", "Seasonal_allergies",
                  "Migraine", "Acid_reflux", "Lung_disease",
                  "Cardiovascular_disease", "Skin_condition", "Obesity")
for (col in disease_cols) {
  if (col %in% names(meta_df)) {
    n <- sum(meta_df[[col]], na.rm = TRUE)
    if (n > 0) {
      cat(sprintf("  %s: %d (%.1f%%)\n", col, n, 100 * n / nrow(meta_df)))
    }
  }
}

## --- Package as a single list ----------------------------------------------

agp_gut <- list(
  counts = counts,
  meta   = meta_df,
  taxa   = taxa,
  source = paste0(
    "dCST-stratified subsample (n=", nrow(counts),
    ") from the American Gut Project (PRJEB11419, AGP-US-2015). ",
    "Original dataset: ~26,900 gut samples, subsampled from a 5,000-sample ",
    "exploratory analysis. Includes all samples in four uncommon demonstration ",
    "dCSTs (Prevotella_7, Pasteurellaceae, Akkermansia, Staphylococcus), then ",
    "fills the remaining slots by a seed-42 simple random sample from the ",
    "eligible background. Phenotypes do not influence selection. ",
    "Eukaryota and Unassigned depth-1 labels are excluded. ",
    "Counts are raw 16S V4 reads (SILVA taxonomy). ",
    "Metadata includes self-reported health conditions from the AGP questionnaire. ",
    "Created with set.seed(42), March 2026. ",
    "Run IDs, selection reasons, and derived annotations are retained in ",
    "inst/extdata/agp_gut_meta.csv. See data-raw/create_agp_gut_subset.py, ",
    "data-raw/build_agp_gut.R, ",
    "and inst/DATA_PROVENANCE.md for provenance and appropriate-use details."
  )
)

## --- Save ------------------------------------------------------------------

save(agp_gut, file = "data/agp_gut.rda", compress = "xz")

cat("\nSaved data/agp_gut.rda\n")
cat("File size:", round(file.info("data/agp_gut.rda")$size / 1024), "KB\n")
