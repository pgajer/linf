#!/usr/bin/env Rscript

# Prepare the pinned public VALENCIA CSVs for the bundled dataset builders.
# Usage: Rscript data-raw/prepare_valencia.R INPUT_DIRECTORY OUTPUT_DIRECTORY
# Input hashes and upstream revision are in inst/DATA_MANIFEST.csv.

normalize_label <- function(x) {
    y <- trimws(as.character(x))
    y[y == ""] <- NA_character_
    y
}

pick_first_non_missing <- function(x) {
    y <- normalize_label(x)
    y <- y[!is.na(y)]
    if (length(y) == 0L) {
        return(NA_character_)
    }
    y[[1L]]
}

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) stop("Usage: prepare_valencia.R INPUT_DIRECTORY OUTPUT_DIRECTORY")
valencia.dir <- normalizePath(args[[1]], mustWork = TRUE)
output.dir <- args[[2]]
dir.create(output.dir, recursive = TRUE, showWarnings = FALSE)
output.dir <- normalizePath(output.dir, mustWork = TRUE)
tax.path <- file.path(valencia.dir, "all_samples_taxonomic_composition_data.csv")
meta.path <- file.path(valencia.dir, "all_samples_metadata.csv")

out.tx.path <- file.path(output.dir, "tx.13k.rds")
out.cst.path <- file.path(output.dir, "cst.tx.13k.rds")
out.summary.path <- file.path(output.dir, "tx.13k_build_summary.rds")
out.sample.filter.csv <- file.path(output.dir, "tx.13k_sample_filtering_summary.csv")
out.taxon.filter.csv <- file.path(output.dir, "tx.13k_taxon_filtering_summary.csv")

MIN_SAMPLE_READS <- 3000L
MIN_COUNT_FOR_PRESENCE <- 2L
MIN_PREVALENCE_PROP <- 0.01
MIN_TOTAL_TAXON_READS <- 50L

if (!file.exists(tax.path)) stop("Missing file: ", tax.path)
if (!file.exists(meta.path)) stop("Missing file: ", meta.path)

manifest <- read.csv("inst/DATA_MANIFEST.csv", stringsAsFactors = FALSE)
raw <- manifest[manifest$kind == "upstream", ]
stopifnot(identical(unname(tools::md5sum(file.path(valencia.dir, raw$file))), raw$md5))

cat("Reading input files...\n")
tax.df <- read.csv(tax.path, check.names = FALSE, stringsAsFactors = FALSE)
meta.df <- read.csv(meta.path, check.names = FALSE, stringsAsFactors = FALSE)

required.tax.cols <- c(
    "Sample_number_for_SRA",
    "Subject_number",
    "HC_CST",
    "HC_subCST",
    "Val_CST",
    "Val_subCST",
    "total_reads"
)
missing.tax.cols <- setdiff(required.tax.cols, colnames(tax.df))
if (length(missing.tax.cols) > 0L) {
    stop("Tax table is missing required columns: ", paste(missing.tax.cols, collapse = ", "))
}

taxa.cols <- setdiff(colnames(tax.df), required.tax.cols)
if (length(taxa.cols) == 0L) {
    stop("No taxa columns found in tax table after removing metadata columns.")
}

tax.counts <- as.matrix(tax.df[, taxa.cols, drop = FALSE])
storage.mode(tax.counts) <- "double"
if (any(!is.finite(tax.counts))) stop("Tax count matrix contains non-finite values.")
if (any(tax.counts < 0)) stop("Tax count matrix contains negative values.")

sample.ids.all <- as.character(tax.df$Sample_number_for_SRA)
rownames(tax.counts) <- sample.ids.all

cat("Filtering samples by minimum read depth before rare-taxa filtering...\n")
sample.reads.pre <- rowSums(tax.counts)
keep.sample.pre <- sample.reads.pre >= MIN_SAMPLE_READS
tax.counts.pre <- tax.counts[keep.sample.pre, , drop = FALSE]

cat("Filtering rare taxa...\n")
taxon.prevalence <- colMeans(tax.counts.pre >= MIN_COUNT_FOR_PRESENCE)
taxon.total.reads <- colSums(tax.counts.pre)
keep.taxa <- (taxon.prevalence >= MIN_PREVALENCE_PROP) &
    (taxon.total.reads >= MIN_TOTAL_TAXON_READS)
tax.counts.filt <- tax.counts.pre[, keep.taxa, drop = FALSE]

cat("Filtering samples by minimum read depth after rare-taxa filtering...\n")
sample.reads.post <- rowSums(tax.counts.filt)
keep.sample.post <- sample.reads.post >= MIN_SAMPLE_READS
tax.counts.final <- tax.counts.filt[keep.sample.post, , drop = FALSE]

if (nrow(tax.counts.final) == 0L) stop("No samples retained after filtering.")
if (ncol(tax.counts.final) == 0L) stop("No taxa retained after filtering.")

row.sums.final <- rowSums(tax.counts.final)
if (any(row.sums.final <= 0)) stop("Found sample(s) with zero total counts after filtering.")

cat("L1-normalizing rows to produce tx.13k...\n")
tx.13k <- tax.counts.final / row.sums.final
storage.mode(tx.13k) <- "double"

## Build metadata-unique lookup by Sample_number.
meta.df$Sample_number <- as.character(meta.df$Sample_number)
meta.df$Subject_number <- as.character(meta.df$Subject_number)
meta.df$CST <- normalize_label(meta.df$CST)
meta.df$subCST <- normalize_label(meta.df$subCST)

meta.idx.split <- split(seq_len(nrow(meta.df)), meta.df$Sample_number)
meta.unique <- do.call(
    rbind,
    lapply(meta.idx.split, function(idx) {
        rows <- meta.df[idx, , drop = FALSE]
        data.frame(
            Sample_number = rows$Sample_number[[1L]],
            Subject_number_meta = rows$Subject_number[[1L]],
            CST_meta = pick_first_non_missing(rows$CST),
            subCST_meta = pick_first_non_missing(rows$subCST),
            stringsAsFactors = FALSE
        )
    })
)
rownames(meta.unique) <- NULL

tax.annotation <- data.frame(
    sample_id = as.character(tax.df$Sample_number_for_SRA),
    Subject_number = as.character(tax.df$Subject_number),
    HC_CST = normalize_label(tax.df$HC_CST),
    HC_subCST = normalize_label(tax.df$HC_subCST),
    Val_CST = normalize_label(tax.df$Val_CST),
    Val_subCST = normalize_label(tax.df$Val_subCST),
    stringsAsFactors = FALSE
)

final.sample.ids <- rownames(tx.13k)
tax.match <- match(final.sample.ids, tax.annotation$sample_id)
if (any(is.na(tax.match))) {
    stop("Internal error: could not align tax annotation to tx.13k samples.")
}
tax.ann.final <- tax.annotation[tax.match, , drop = FALSE]

meta.match <- match(final.sample.ids, meta.unique$Sample_number)
meta.cst <- meta.unique$CST_meta[meta.match]
meta.subcst <- meta.unique$subCST_meta[meta.match]

cst.final <- ifelse(!is.na(meta.cst), meta.cst, tax.ann.final$Val_CST)
subcst.final <- ifelse(!is.na(meta.subcst), meta.subcst, tax.ann.final$Val_subCST)

cst.tx.13k <- data.frame(
    sample_id = final.sample.ids,
    Subject_number = tax.ann.final$Subject_number,
    CST = cst.final,
    subCST = subcst.final,
    CST_source = ifelse(!is.na(meta.cst), "metadata_CST", "Val_CST"),
    subCST_source = ifelse(!is.na(meta.subcst), "metadata_subCST", "Val_subCST"),
    metadata_CST = meta.cst,
    metadata_subCST = meta.subcst,
    Val_CST = tax.ann.final$Val_CST,
    Val_subCST = tax.ann.final$Val_subCST,
    HC_CST = tax.ann.final$HC_CST,
    HC_subCST = tax.ann.final$HC_subCST,
    stringsAsFactors = FALSE
)
rownames(cst.tx.13k) <- cst.tx.13k$sample_id

sample.filter.summary <- data.frame(
    sample_id = sample.ids.all,
    reads_pre_filter = sample.reads.pre,
    keep_pre_sample_filter = keep.sample.pre,
    stringsAsFactors = FALSE
)
sample.filter.summary <- sample.filter.summary[match(rownames(tax.counts.pre), sample.filter.summary$sample_id), , drop = FALSE]
sample.filter.summary$reads_after_taxa_filter <- rowSums(tax.counts.filt)
sample.filter.summary$keep_post_sample_filter <- keep.sample.post
sample.filter.summary$kept_in_tx13k <- sample.filter.summary$sample_id %in% final.sample.ids

taxon.filter.summary <- data.frame(
    taxon = colnames(tax.counts.pre),
    prevalence = as.numeric(taxon.prevalence),
    total_reads = as.numeric(taxon.total.reads),
    keep_taxon = as.logical(keep.taxa),
    stringsAsFactors = FALSE
)

summary.obj <- list(
    params = list(
        min.sample.reads = MIN_SAMPLE_READS,
        min.count.for.presence = MIN_COUNT_FOR_PRESENCE,
        min.prevalence.prop = MIN_PREVALENCE_PROP,
        min.total.taxon.reads = MIN_TOTAL_TAXON_READS
    ),
    dimensions = list(
        n.samples.input = nrow(tax.counts),
        n.taxa.input = ncol(tax.counts),
        n.samples.after.pre.sample.filter = nrow(tax.counts.pre),
        n.taxa.after.taxa.filter = ncol(tax.counts.filt),
        n.samples.final = nrow(tx.13k),
        n.taxa.final = ncol(tx.13k)
    ),
    missing.labels = list(
        cst_missing = sum(is.na(cst.tx.13k$CST)),
        subcst_missing = sum(is.na(cst.tx.13k$subCST))
    )
)

saveRDS(tx.13k, out.tx.path, compress = "xz")
saveRDS(cst.tx.13k, out.cst.path, compress = "xz")
saveRDS(summary.obj, out.summary.path, compress = "xz")
write.csv(sample.filter.summary, out.sample.filter.csv, row.names = FALSE)
write.csv(taxon.filter.summary, out.taxon.filter.csv, row.names = FALSE)

cat("\nDone.\n")
cat("tx.13k samples:", nrow(tx.13k), " taxa:", ncol(tx.13k), "\n")
cat("Saved tx.13k:", out.tx.path, "\n")
cat("Saved CST table:", out.cst.path, "\n")
cat("Saved summary:", out.summary.path, "\n")
