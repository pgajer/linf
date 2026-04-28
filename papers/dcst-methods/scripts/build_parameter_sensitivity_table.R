#!/usr/bin/env Rscript

`%||%` <- function(x, y) if (is.null(x)) y else x

script_path <- sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1] %||% "")
script_dir <- if (nzchar(script_path)) dirname(normalizePath(script_path, mustWork = TRUE)) else getwd()
repo <- normalizePath(file.path(script_dir, "../../.."), mustWork = FALSE)
if (!dir.exists(file.path(repo, "papers", "dcst-methods"))) {
  repo <- normalizePath(getwd(), mustWork = TRUE)
}

fmt_pct <- function(x, digits = 1) sprintf(paste0("%.", digits, "f\\%%"), 100 * x)
fmt_q <- function(x) {
  ifelse(is.na(x), "NA", ifelse(x < 0.001, sprintf("%.2e", x), sprintf("%.4f", x)))
}
fmt_v <- function(x) sprintf("%.3f", x)

tables_dir <- file.path(repo, "papers", "dcst-methods", "assets", "tables")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

raw_sizes_path <- file.path(
  repo,
  "papers", "gut-dcst", "notes", "n0_quantile_threshold_report",
  "tables", "agp_raw_lineage_set_sizes_by_depth.tsv"
)
ct_cap_path <- file.path(
  "/Users/pgajer/current_projects/CT_clearance",
  "analysis_output/ct_combined_vog_dcst_sensitivity_2026-04-10",
  "sensitivity_cap_summary.tsv"
)
agp_demo_path <- file.path(repo, "data", "agp_gut.rda")
ct_cap50_matrix_path <- file.path(
  "/Users/pgajer/current_projects/CT_clearance",
  "data/processed/ct_biological_plus_virgo2_sensitivity_2026-04-10",
  "cap_050000/x__ct_biological_plus_virgo2__vog.rds"
)

stopifnot(file.exists(raw_sizes_path), file.exists(ct_cap_path), file.exists(agp_demo_path))

raw_sizes <- read.delim(raw_sizes_path, stringsAsFactors = FALSE)
raw_sizes <- raw_sizes[raw_sizes$depth %in% 1:4, , drop = FALSE]
nominal <- c(`1` = 50L, `2` = 25L, `3` = 25L, `4` = 25L)
variants <- data.frame(
  variant = c("-25\\%", "nominal", "+25\\%"),
  factor = c(0.75, 1.00, 1.25),
  stringsAsFactors = FALSE
)
support_rows <- do.call(rbind, lapply(seq_len(nrow(variants)), function(i) {
  do.call(rbind, lapply(1:4, function(d) {
    dat <- raw_sizes[raw_sizes$depth == d, , drop = FALSE]
    threshold <- floor(nominal[as.character(d)] * variants$factor[i] + 0.5)
    keep <- dat$size >= threshold
    data.frame(
      variant = variants$variant[i],
      depth = d,
      threshold = threshold,
      retained_lineages = sum(keep),
      retained_share = sum(dat$size[keep]) / sum(dat$size),
      stringsAsFactors = FALSE
    )
  }))
}))
support_tsv <- file.path(tables_dir, "TABLE_S4A_agp_support_schedule_perturbation.tsv")
write.table(support_rows, support_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

support_lines <- vapply(split(support_rows, support_rows$variant), function(dat) {
  dat <- dat[order(dat$depth), ]
  paste0(
    dat$variant[1], ": thresholds ",
    paste(dat$threshold, collapse = "/"),
    "; retained lineages ",
    paste(dat$retained_lineages, collapse = "/"),
    "; pre-absorb retained sample share ",
    paste(fmt_pct(dat$retained_share), collapse = "/")
  )
}, character(1))
support_summary <- paste(support_lines[c("-25\\%", "nominal", "+25\\%")], collapse = ". ")

tie_audit <- function(X) {
  mx <- apply(X, 1, max)
  ties <- rowSums(X == mx)
  data.frame(
    n_samples = nrow(X),
    n_features = ncol(X),
    zero_rows = sum(mx <= 0),
    top_tie_rows = sum(mx > 0 & ties > 1),
    top_tie_pct = mean(mx > 0 & ties > 1),
    max_top_tie_size = max(ties[mx > 0]),
    stringsAsFactors = FALSE
  )
}

load(agp_demo_path)
agp_tie <- tie_audit(agp_gut$counts)
ct_tie_note <- "CT tie audit not rerun"
ct_tie <- NULL
if (file.exists(ct_cap50_matrix_path)) {
  ct_cap50 <- readRDS(ct_cap50_matrix_path)
  ct_tie <- tie_audit(ct_cap50)
  rm(ct_cap50)
  gc()
  ct_tie_note <- paste0(
    ct_tie$top_tie_rows, "/", ct_tie$n_samples,
    " rows in the 50k-VOG cap audit (", fmt_pct(ct_tie$top_tie_pct, 2),
    "); maximum tied top set size ", ct_tie$max_top_tie_size
  )
}
tie_rows <- rbind(
  data.frame(dataset = "Bundled AGP gut example", agp_tie, stringsAsFactors = FALSE),
  if (!is.null(ct_tie)) data.frame(dataset = "CT 50k-VOG cap matrix", ct_tie, stringsAsFactors = FALSE)
)
tie_tsv <- file.path(tables_dir, "TABLE_S4B_top_dominance_tie_audit.tsv")
write.table(tie_rows, tie_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

ct_cap <- read.delim(ct_cap_path, stringsAsFactors = FALSE)
ct_cap <- ct_cap[ct_cap$cap %in% c(25000, 50000, 100000, 250000, 509501), ]
ct_cap <- ct_cap[order(ct_cap$cap), ]
ct_summary <- paste0(
  paste0(
    ct_cap$cap_label,
    ": depth-3 q=", fmt_q(ct_cap$depth3_omnibus_q),
    ", V=", fmt_v(ct_cap$depth3_cramers_v),
    ", depth-4 q=", fmt_q(ct_cap$depth4_omnibus_q)
  ),
  collapse = "; "
)

rows <- data.frame(
  Check = c(
    "Deterministic tie rule",
    "Observed top-feature ties",
    "Gut support-schedule perturbation",
    "CT retained-feature cap robustness",
    "Gut feature-filter scope"
  ),
  Summary = c(
    paste0(
      "All reported dCST fits used \\texttt{tie.method = \"first\"}: first retained feature for primary dominance and first retained target for absorb reassignment. ",
      "This is part of the label definition, not a random analysis option."
    ),
    paste0(
      "Dominant-feature ties were absent in the bundled AGP gut example (",
      agp_tie$top_tie_rows, "/", agp_tie$n_samples,
      " rows) and rare in the CT cap audit: ", ct_tie_note,
      "."
    ),
    paste0(
      "Using the existing full-AGP raw-lineage size audit, a \\(\\pm25\\%\\) perturbation around the manuscript support schedule 50/25/25/25 changed hierarchy granularity but preserved high pre-absorb sample coverage at shallow depths. ",
      support_summary,
      "."
    ),
    paste0(
      "Capped VOG analyses used nested variance-ranked subsets from the same 509,501-VOG retained universe. ",
      ct_summary,
      ". The 25k cap was shallower, whereas 50k and larger caps retained depth-4 support."
    ),
    paste0(
      "The gut taxonomic analysis is conditional on the stated shared-taxonomy retained family: count \\(\\geq2\\) in at least 5\\% of retained AGP samples, yielding 340 taxa. ",
      "The combined manuscript does not claim robustness to all alternative gut prevalence filters."
    )
  ),
  Interpretation = c(
    "The reported labels are reproducible by definition once feature order is fixed.",
    "Tie handling is specified for exact reproducibility; available audits do not suggest that top-feature ties drive the reported examples.",
    "The support schedule controls naming granularity and should be treated as an explicit modeling choice, not as a hidden clustering parameter.",
    "This supports the CT result as a within-cohort feature-family robustness result, not as external replication.",
    "This bounds the gut portability claim to the retained feature family used for transfer and rebuilt analyses."
  ),
  stringsAsFactors = FALSE
)

summary_tsv <- file.path(tables_dir, "TABLE_S4_parameter_sensitivity_summary.tsv")
write.table(rows, summary_tsv, sep = "\t", quote = FALSE, row.names = FALSE)

escape_latex <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x
}

tex_path <- file.path(tables_dir, "TABLE_S4_parameter_sensitivity_summary.tex")
con <- file(tex_path, open = "wt")
on.exit(close(con), add = TRUE)
writeLines(c(
  "\\clearpage",
  "",
  "\\noindent\\textbf{Supplementary Table S4.} Setting-dependence and robustness",
  "checks for deterministic dCST construction.",
  "",
  "\\begin{longtable}{p{0.18\\linewidth}p{0.48\\linewidth}p{0.26\\linewidth}}",
  "\\toprule",
  "Check & Summary & Interpretation \\\\",
  "\\midrule",
  "\\endfirsthead",
  "\\toprule",
  "Check & Summary & Interpretation \\\\",
  "\\midrule",
  "\\endhead",
  "\\bottomrule",
  "\\endfoot"
), con)
for (i in seq_len(nrow(rows))) {
  writeLines(paste0(
    escape_latex(rows$Check[i]), " & ",
    rows$Summary[i], " & ",
    rows$Interpretation[i], " \\\\"
  ), con)
  if (i < nrow(rows)) writeLines("\\midrule", con)
}
writeLines("\\end{longtable}", con)

message("Wrote ", tex_path)
