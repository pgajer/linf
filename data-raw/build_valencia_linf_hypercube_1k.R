## ============================================================================
## build_valencia_linf_hypercube_1k.R
## ============================================================================
##
## Creates the bundled `valencia_linf_hypercube_1k` example dataset from the
## canonical Valencia 13k training set. The object is intentionally small and
## focused: it contains four biologically interpretable phylotype coordinates
## for demonstrating the zero-aware homogeneous-coordinate hypercube embedding.
##
## Input files (not shipped with the package):
##   $LINF_VALENCIA_ROOT/tx.13k.rds      -- compositional matrix
##   $LINF_VALENCIA_ROOT/cst.tx.13k.rds  -- CST annotation data frame
##
## Output (saved to data/):
##   valencia_linf_hypercube_1k
##
## Run from the package root:
##   source("data-raw/build_valencia_linf_hypercube_1k.R")
## ============================================================================

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

component_map <- c(
  Li = "Lactobacillus_iners",
  Lc = "Lactobacillus_crispatus",
  Gv = "Gardnerella_vaginalis",
  Bv = "BVAB1"
)

missing_taxa <- setdiff(unname(component_map), colnames(tx))
if (length(missing_taxa)) {
  stop("Missing required Valencia taxa: ", paste(missing_taxa, collapse = ", "))
}

rel4_raw <- as.matrix(tx[, unname(component_map), drop = FALSE])
colnames(rel4_raw) <- names(component_map)

selected_mass <- rowSums(rel4_raw)
keep <- selected_mass > 0
rel4_raw <- rel4_raw[keep, , drop = FALSE]
cst_keep <- cst_df[keep, , drop = FALSE]
source_rows <- which(keep)
selected_mass <- selected_mass[keep]

rel4 <- sweep(rel4_raw, 1L, selected_mass, "/")
dominant_component <- colnames(rel4)[max.col(rel4, ties.method = "first")]

allocate_proportional <- function(tab, n_target) {
  raw <- as.numeric(tab) / sum(tab) * n_target
  alloc <- floor(raw)
  shortfall <- n_target - sum(alloc)
  if (shortfall > 0L) {
    ord <- order(raw - alloc, decreasing = TRUE)
    alloc[ord[seq_len(shortfall)]] <- alloc[ord[seq_len(shortfall)]] + 1L
  }
  if (any(alloc > as.numeric(tab))) {
    stop("Requested allocation exceeds available rows in a stratum")
  }
  names(alloc) <- names(tab)
  alloc
}

set.seed(20261604L)
n_target <- 1000L
tab <- table(dominant_component)
n_by_group <- allocate_proportional(tab, n_target)

sampled <- integer(0)
for (group in names(n_by_group)) {
  idx <- which(dominant_component == group)
  sampled <- c(sampled, sample(idx, n_by_group[[group]]))
}
sampled <- sort(sampled)

rel4 <- rel4[sampled, , drop = FALSE]
cst_sub <- cst_keep[sampled, , drop = FALSE]
source_rows <- source_rows[sampled]
selected_mass <- selected_mass[sampled]
dominant_component <- dominant_component[sampled]

sample_id <- sprintf("v1k_%04d", seq_len(nrow(rel4)))
rownames(rel4) <- sample_id

meta <- data.frame(
  sample_id = sample_id,
  source_row = source_rows,
  Val_CST = as.character(cst_sub$Val_CST),
  Val_subCST = as.character(cst_sub$Val_subCST),
  selected_mass = as.numeric(selected_mass),
  dominant_component = factor(dominant_component, levels = names(component_map)),
  stringsAsFactors = FALSE
)

stopifnot(nrow(rel4) == n_target)
stopifnot(all(abs(rowSums(rel4) - 1) < 1e-10))
stopifnot(identical(rownames(rel4), meta$sample_id))

cat("Dominant component distribution:\n")
print(table(meta$dominant_component))

cat("\nValencia CST distribution:\n")
print(sort(table(meta$Val_CST), decreasing = TRUE))

valencia_linf_hypercube_1k <- list(
  rel4 = rel4,
  meta = meta,
  component_map = component_map,
  source = paste0(
    "Stratified n=1000 four-component subset from the Valencia 13k vaginal ",
    "16S training set. Components are Li=Lactobacillus_iners, ",
    "Lc=Lactobacillus_crispatus, Gv=Gardnerella_vaginalis, and Bv=BVAB1. ",
    "Rows with zero mass in these four components were removed, remaining ",
    "rows were L1-normalized over the selected components, then sampled ",
    "with set.seed(20261604), stratified by dominant selected component."
  )
)

save(
  valencia_linf_hypercube_1k,
  file = "data/valencia_linf_hypercube_1k.rda",
  compress = "xz"
)

cat("\nSaved data/valencia_linf_hypercube_1k.rda\n")
cat("File size:", round(file.info("data/valencia_linf_hypercube_1k.rda")$size / 1024), "KB\n")
