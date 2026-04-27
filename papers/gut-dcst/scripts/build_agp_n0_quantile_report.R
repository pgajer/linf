#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(scales)
})

stopf <- function(...) stop(sprintf(...), call. = FALSE)
messagef <- function(...) message(sprintf(...))

repo_root <- "/Users/pgajer/current_projects/linf"
gut_root <- "/Users/pgajer/current_projects/gut_microbiome"
counts_path <- file.path(
  gut_root,
  "outputs/prime_species/prime_gut_projects_silva_species_absolute_2026-03-24.csv.gz"
)
metadata_path <- file.path(
  gut_root,
  "data/prime_gut_project_sample_metadata_2026-03-24.csv.gz"
)
reference_run_dir <- file.path(
  gut_root,
  "outputs/dcst_analysis/runs/2026-04-11-absorb-depthscan-adaptive"
)
out_dir <- file.path(repo_root, "papers/gut-dcst/notes/n0_quantile_threshold_report")
fig_dir <- file.path(out_dir, "figures")
table_dir <- file.path(out_dir, "tables")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

current_n0 <- c(50L, 25L, 12L, 10L, 8L, 5L)
candidate_q <- c(0.50, 0.60, 0.70, 0.75, 0.80, 0.90, 0.95, 0.975, 0.99)
max_depth <- 6L

load_linf <- function(linf_root) {
  if (requireNamespace("pkgload", quietly = TRUE)) {
    pkgload::load_all(linf_root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
    return(invisible(TRUE))
  }
  if (requireNamespace("devtools", quietly = TRUE)) {
    devtools::load_all(linf_root, quiet = TRUE, helpers = FALSE, export_all = FALSE)
    return(invisible(TRUE))
  }
  stopf("Neither pkgload nor devtools is available; cannot load %s", linf_root)
}

read_metadata <- function(path, agp_project = "PRJEB11419") {
  cols <- c("Run", "BioProject")
  meta <- fread(cmd = sprintf("gzip -dc %s", shQuote(path)), select = cols, data.table = FALSE)
  meta <- meta[meta$BioProject == agp_project, , drop = FALSE]
  if (!nrow(meta)) stopf("No metadata rows found for project %s", agp_project)
  meta
}

read_counts_table <- function(path, runs_to_keep) {
  counts_dt <- fread(cmd = sprintf("gzip -dc %s", shQuote(path)))
  counts_dt <- counts_dt[Run %in% runs_to_keep]
  if (!nrow(counts_dt)) stopf("No count rows matched the requested AGP runs")
  run_ids <- counts_dt$Run
  count_mat <- as.matrix(counts_dt[, setdiff(names(counts_dt), "Run"), with = FALSE])
  storage.mode(count_mat) <- "integer"
  rownames(count_mat) <- run_ids
  rm(counts_dt)
  invisible(gc())
  count_mat
}

extract_rank_value <- function(parts, prefix) {
  hit <- parts[startsWith(parts, prefix)]
  if (!length(hit)) return("")
  value <- sub(paste0("^", prefix), "", hit[[1L]])
  value <- trimws(value)
  if (!nzchar(value) || grepl("^_+$", value)) return("")
  value
}

compact_taxon_label <- function(taxonomy) {
  if (is.na(taxonomy) || !nzchar(taxonomy)) return("Unassigned")
  parts <- strsplit(taxonomy, ";", fixed = TRUE)[[1L]]
  species <- extract_rank_value(parts, "s__")
  genus <- extract_rank_value(parts, "g__")
  family <- extract_rank_value(parts, "f__")
  order <- extract_rank_value(parts, "o__")
  class_name <- extract_rank_value(parts, "c__")
  phylum <- extract_rank_value(parts, "p__")
  domain <- extract_rank_value(parts, "d__")

  if (nzchar(species)) {
    if (nzchar(genus) && !startsWith(species, genus)) return(paste(genus, species))
    return(species)
  }
  if (nzchar(genus)) return(genus)
  if (nzchar(family)) return(family)
  if (nzchar(order)) return(order)
  if (nzchar(class_name)) return(class_name)
  if (nzchar(phylum)) return(phylum)
  if (nzchar(domain)) return(domain)
  "Unassigned"
}

build_feature_reference <- function(raw_taxa) {
  base_labels <- vapply(raw_taxa, compact_taxon_label, character(1L))
  data.frame(
    feature_id = sprintf("tx%04d", seq_along(raw_taxa)),
    feature_label = make.unique(base_labels, sep = "_"),
    taxonomy = raw_taxa,
    stringsAsFactors = FALSE
  )
}

build_rank_prefix_matrix <- function(counts, feature_labels, max_depth = 6L) {
  n <- nrow(counts)
  prefixes <- matrix(NA_character_, nrow = n, ncol = max_depth)
  rownames(prefixes) <- rownames(counts)
  colnames(prefixes) <- sprintf("depth%d", seq_len(max_depth))

  for (i in seq_len(n)) {
    vals <- counts[i, ]
    pos <- which(vals > 0)
    if (!length(pos)) next
    ord <- pos[order(-vals[pos], pos)]
    keep <- ord[seq_len(min(max_depth, length(ord)))]
    labs <- feature_labels[keep]
    prefixes[i, seq_along(labs)] <- vapply(seq_along(labs), function(j) {
      paste(labs[seq_len(j)], collapse = "__")
    }, character(1L))
  }

  prefixes
}

summarize_size_distributions <- function(prefixes, current_n0) {
  pieces <- vector("list", ncol(prefixes))
  for (d in seq_len(ncol(prefixes))) {
    lineage <- prefixes[, d]
    lineage <- lineage[!is.na(lineage)]
    tab <- sort(table(lineage), decreasing = TRUE)
    sizes <- as.integer(tab)
    qs <- quantile(sizes, probs = c(0.10, 0.25, 0.50, 0.75, 0.80, 0.90, 0.95, 0.99), type = 7)
    pieces[[d]] <- data.frame(
      depth = d,
      n_samples_with_depth = length(lineage),
      n_lineage_sets = length(sizes),
      min_size = min(sizes),
      p10 = unname(qs[["10%"]]),
      p25 = unname(qs[["25%"]]),
      median = unname(qs[["50%"]]),
      mean = mean(sizes),
      p75 = unname(qs[["75%"]]),
      p80 = unname(qs[["80%"]]),
      p90 = unname(qs[["90%"]]),
      p95 = unname(qs[["95%"]]),
      p99 = unname(qs[["99%"]]),
      max_size = max(sizes),
      current_n0 = current_n0[[d]],
      current_ecdf = mean(sizes <= current_n0[[d]]),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pieces)
}

candidate_thresholds <- function(prefixes, probs) {
  pieces <- list()
  k <- 1L
  for (d in seq_len(ncol(prefixes))) {
    lineage <- prefixes[, d]
    lineage <- lineage[!is.na(lineage)]
    sizes <- as.integer(table(lineage))
    for (q in probs) {
      n0 <- ceiling(unname(quantile(sizes, probs = q, type = 7)))
      keep <- sizes >= n0
      pieces[[k]] <- data.frame(
        q = q,
        depth = d,
        n0 = n0,
        n_retained_lineage_sets = sum(keep),
        retained_sample_share = sum(sizes[keep]) / sum(sizes),
        stringsAsFactors = FALSE
      )
      k <- k + 1L
    }
  }
  do.call(rbind, pieces)
}

lineage_size_long <- function(prefixes) {
  pieces <- vector("list", ncol(prefixes))
  for (d in seq_len(ncol(prefixes))) {
    lineage <- prefixes[, d]
    lineage <- lineage[!is.na(lineage)]
    tab <- table(lineage)
    pieces[[d]] <- data.frame(
      depth = d,
      lineage = names(tab),
      size = as.integer(tab),
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, pieces)
}

latex_escape <- function(x) {
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([#$%&_{}])", "\\\\\\1", x, perl = TRUE)
  x <- gsub("~", "\\\\textasciitilde{}", x, fixed = TRUE)
  x <- gsub("\\^", "\\\\textasciicircum{}", x)
  x
}

format_num <- function(x, digits = 1) {
  ifelse(abs(x - round(x)) < 1e-8, as.character(round(x)), sprintf(paste0("%.", digits, "f"), x))
}

write_latex_table <- function(df, path, align = NULL) {
  if (is.null(align)) align <- paste0("l", paste(rep("r", ncol(df) - 1L), collapse = ""))
  con <- file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  writeLines(sprintf("\\begin{tabular}{%s}", align), con)
  writeLines("\\toprule", con)
  writeLines(paste(latex_escape(names(df)), collapse = " & "), con)
  writeLines(" \\\\", con)
  writeLines("\\midrule", con)
  for (i in seq_len(nrow(df))) {
    row <- vapply(df[i, ], as.character, character(1L))
    writeLines(paste(latex_escape(row), collapse = " & "), con)
    writeLines(" \\\\", con)
  }
  writeLines("\\bottomrule", con)
  writeLines("\\end{tabular}", con)
}

message("Loading linf package functions...")
load_linf(repo_root)

message("Loading AGP metadata and PRIME SILVA count table...")
meta <- read_metadata(metadata_path)
counts_matrix <- read_counts_table(counts_path, runs_to_keep = meta$Run)
raw_n_samples <- nrow(counts_matrix)
raw_n_taxa <- ncol(counts_matrix)

message("Applying the discovery filters used in the absorb depth-scan...")
filt <- filter.asv(
  counts_matrix,
  min.lib = 1000,
  prev.prop = 0.05,
  min.count = 2,
  verbose = TRUE
)
rm(counts_matrix)
invisible(gc())

feature_reference <- build_feature_reference(colnames(filt$counts))
colnames(filt$counts) <- feature_reference$feature_label

message("Reconstructing pre-absorption dominance-lineage rank prefixes...")
prefixes <- build_rank_prefix_matrix(filt$counts, feature_reference$feature_label, max_depth = max_depth)
lineage_sizes <- lineage_size_long(prefixes)
size_summary <- summarize_size_distributions(prefixes, current_n0 = current_n0)
thresholds <- candidate_thresholds(prefixes, candidate_q)

message("Writing analysis tables...")
fwrite(size_summary, file.path(table_dir, "agp_raw_lineage_size_summary.tsv"), sep = "\t")
fwrite(thresholds, file.path(table_dir, "agp_quantile_n0_candidates.tsv"), sep = "\t")
fwrite(lineage_sizes[order(lineage_sizes$depth, -lineage_sizes$size), ],
       file.path(table_dir, "agp_raw_lineage_set_sizes_by_depth.tsv"), sep = "\t")
fwrite(feature_reference, file.path(table_dir, "agp_filtered_feature_reference.tsv"), sep = "\t")

threshold_wide <- dcast(as.data.table(thresholds), q ~ depth, value.var = "n0")
setnames(threshold_wide, old = as.character(seq_len(max_depth)), new = paste0("depth_", seq_len(max_depth)))
threshold_wide[, q := sprintf("%.3f", as.numeric(q))]
fwrite(threshold_wide, file.path(table_dir, "agp_quantile_n0_candidates_wide.tsv"), sep = "\t")

retention_wide <- as.data.table(thresholds)
retention_wide[, retained_sample_share := sprintf("%.1f%%", 100 * retained_sample_share)]
retention_wide <- dcast(as.data.table(retention_wide), q ~ depth, value.var = "retained_sample_share")
setnames(retention_wide, old = as.character(seq_len(max_depth)), new = paste0("depth_", seq_len(max_depth)))
retention_wide[, q := sprintf("%.3f", as.numeric(q))]
fwrite(retention_wide, file.path(table_dir, "agp_quantile_retained_sample_share_wide.tsv"), sep = "\t")

summary_tex <- data.frame(
  Depth = size_summary$depth,
  `Lineage sets` = size_summary$n_lineage_sets,
  Min = size_summary$min_size,
  Median = format_num(size_summary$median),
  Mean = format_num(size_summary$mean),
  P75 = format_num(size_summary$p75),
  P90 = format_num(size_summary$p90),
  P95 = format_num(size_summary$p95),
  Max = size_summary$max_size,
  `Current n0` = size_summary$current_n0,
  `Current ECDF` = sprintf("%.1f%%", 100 * size_summary$current_ecdf),
  check.names = FALSE
)
write_latex_table(summary_tex, file.path(table_dir, "lineage_size_summary.tex"))
write_latex_table(as.data.frame(threshold_wide), file.path(table_dir, "quantile_n0_candidates.tex"))
write_latex_table(as.data.frame(retention_wide), file.path(table_dir, "quantile_retained_sample_share.tex"))

message("Writing histograms...")
vline_df <- merge(
  data.frame(depth = seq_len(max_depth), current_n0 = current_n0),
  thresholds[thresholds$q %in% c(0.75, 0.90), c("q", "depth", "n0")],
  by = "depth",
  all.x = TRUE
)

combined_hist <- ggplot(lineage_sizes, aes(x = size)) +
  geom_histogram(bins = 45, fill = "#4E79A7", color = "white", linewidth = 0.2) +
  geom_vline(
    data = data.frame(depth = seq_len(max_depth), value = current_n0, line = "Current n0"),
    aes(xintercept = value, linetype = line),
    color = "#E15759",
    linewidth = 0.55
  ) +
  geom_vline(
    data = thresholds[thresholds$q %in% c(0.75, 0.90), ],
    aes(xintercept = n0, linetype = sprintf("q = %.2f", q)),
    color = "#59A14F",
    linewidth = 0.45
  ) +
  scale_x_log10(labels = label_number()) +
  scale_y_continuous(labels = label_number()) +
  scale_linetype_manual(values = c("Current n0" = "solid", "q = 0.75" = "dashed", "q = 0.90" = "dotdash")) +
  facet_wrap(~ depth, scales = "free_y", ncol = 3, labeller = label_both) +
  labs(
    x = "Lineage-set size (samples, log10 scale)",
    y = "Number of lineage sets",
    linetype = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    strip.background = element_rect(fill = "#F2F2F2", color = NA)
  )
ggsave(file.path(fig_dir, "agp_lineage_size_histograms_depth1_to6.png"),
       combined_hist, width = 8.2, height = 5.8, dpi = 320)

for (d in seq_len(max_depth)) {
  dd <- lineage_sizes[lineage_sizes$depth == d, , drop = FALSE]
  vd <- rbind(
    data.frame(value = current_n0[[d]], label = "Current n0", stringsAsFactors = FALSE),
    data.frame(
      value = thresholds$n0[thresholds$depth == d & thresholds$q %in% c(0.75, 0.90)],
      label = sprintf("q = %.2f", thresholds$q[thresholds$depth == d & thresholds$q %in% c(0.75, 0.90)]),
      stringsAsFactors = FALSE
    )
  )
  p <- ggplot(dd, aes(x = size)) +
    geom_histogram(bins = 45, fill = "#4E79A7", color = "white", linewidth = 0.2) +
    geom_vline(data = vd, aes(xintercept = value, linetype = label),
               color = "#E15759", linewidth = 0.5) +
    scale_x_log10(labels = label_number()) +
    scale_y_continuous(labels = label_number()) +
    labs(
      title = sprintf("Depth %d pre-absorption lineage-set sizes", d),
      x = "Lineage-set size (samples, log10 scale)",
      y = "Number of lineage sets",
      linetype = NULL
    ) +
    theme_bw(base_size = 11) +
    theme(legend.position = "bottom", panel.grid.minor = element_blank())
  ggsave(file.path(fig_dir, sprintf("agp_lineage_size_histogram_depth%d.png", d)),
         p, width = 6.2, height = 4.1, dpi = 320)
}

message("Writing LaTeX report...")
report_tex <- file.path(out_dir, "agp_n0_quantile_threshold_report.tex")
report_lines <- c(
  "\\documentclass[11pt]{article}",
  "\\usepackage[margin=1in]{geometry}",
  "\\usepackage{graphicx}",
  "\\usepackage{booktabs}",
  "\\usepackage{array}",
  "\\usepackage{float}",
  "\\usepackage{placeins}",
  "\\usepackage[T1]{fontenc}",
  "\\usepackage{lmodern}",
  "\\usepackage[hidelinks]{hyperref}",
  "\\title{Data-adaptive minimum-size thresholds for AGP absorb-policy dCST construction}",
  "\\author{Generated from PRIME AGP SILVA discovery data}",
  sprintf("\\date{%s}", format(Sys.Date(), "%B %d, %Y")),
  "\\begin{document}",
  "\\maketitle",
  "",
  "\\section*{Methods-oriented summary}",
  sprintf(
    paste0(
      "This report evaluates data-adaptive candidates for the depth-specific minimum ",
      "lineage-set size parameter $n_0$ used in absorb-policy dCST construction. ",
      "The analysis uses the same AGP discovery source as the current manuscript: ",
      "%s AGP metadata rows were read from the PRIME gut metadata table, the raw ",
      "PRIME SILVA species table contained %s samples by %s taxa before filtering, ",
      "and the discovery filters retained %s samples and %s taxa."
    ),
    format(nrow(meta), big.mark = ","),
    format(raw_n_samples, big.mark = ","),
    format(raw_n_taxa, big.mark = ","),
    format(nrow(filt$counts), big.mark = ","),
    format(ncol(filt$counts), big.mark = ",")
  ),
  "",
  "The purpose is not to replace the planned sensitivity analysis, but to make the default threshold schedule defensible and reproducible. The current manuscript uses $n_0 = 50, 25, 12, 10, 8, 5$ at depths 1 through 6. Here, each candidate schedule is defined by a single quantile $q$: at each depth, $n_0(d)$ is set to the ceiling of the empirical $q$-quantile of the pre-absorption lineage-set size distribution at that depth.",
  "",
  "\\section*{Lineage-set reconstruction}",
  "After applying the established discovery filters (minimum library size 1000 reads; feature prevalence defined as count $\\geq 2$ in at least 5\\% of retained samples), taxa were compacted to the same display labels used by the absorb depth-scan. For each sample, retained taxa were ranked by abundance with deterministic column-order tie handling. A depth-$d$ raw dominance lineage was the ordered prefix of the top $d$ positive-abundance retained taxa. The histograms and tables below summarize these raw rank-prefix sets before low-frequency absorption. This avoids choosing new thresholds from labels already shaped by the current ad hoc $n_0$ schedule.",
  "",
  "\\section*{Empirical lineage-set size distributions}",
  "\\begin{table}[H]",
  "\\centering",
  "\\scriptsize",
  "\\caption{Pre-absorption AGP SILVA dominance-lineage set size distributions by depth. Current ECDF is the empirical percentile position of the current manuscript threshold, computed as the fraction of raw lineage sets with size less than or equal to the current $n_0$.}",
  "\\resizebox{\\linewidth}{!}{\\input{tables/lineage_size_summary.tex}}",
  "\\end{table}",
  "",
  "\\begin{figure}[H]",
  "\\centering",
  "\\includegraphics[width=0.98\\linewidth]{figures/agp_lineage_size_histograms_depth1_to6.png}",
  "\\caption{Histograms of pre-absorption dominance-lineage set sizes for depths 1 through 6. The x-axis is log-scaled because the distributions are strongly right-skewed. Red vertical lines show the current manuscript $n_0$ values; green lines show the 0.75 and 0.90 quantile-derived alternatives.}",
  "\\end{figure}",
  "",
  "\\FloatBarrier",
  "\\section*{Quantile-derived candidate schedules}",
  "Table 2 gives candidate $n_0$ schedules from several plausible quantile choices. Lower values of $q$ preserve more named lineages and push less mass into absorption, whereas higher values of $q$ produce a smaller named hierarchy and more aggressive tail absorption. Table 3 gives the corresponding pre-absorption sample share that would lie in lineage sets at or above the candidate threshold, before applying absorb reassignment. Because raw high-depth rank-prefixes are very sparse, high quantiles are included in addition to more conventional central quantiles.",
  "",
  "\\begin{table}[H]",
  "\\centering",
  "\\scriptsize",
  "\\caption{Candidate $n_0$ values obtained by setting depth-specific $n_0$ to the ceiling of the empirical lineage-set size quantile.}",
  "\\resizebox{\\linewidth}{!}{\\input{tables/quantile_n0_candidates.tex}}",
  "\\end{table}",
  "",
  "\\begin{table}[H]",
  "\\centering",
  "\\scriptsize",
  "\\caption{Pre-absorption sample share in lineage sets retained by each quantile-derived threshold.}",
  "\\resizebox{\\linewidth}{!}{\\input{tables/quantile_retained_sample_share.tex}}",
  "\\end{table}",
  "",
  "\\section*{Interpretation for a methods paper}",
  "The empirical distributions are heavy-tailed at every depth: many lineage sets are small, but a small number of dominant rank-prefixes account for a large fraction of samples. The current manuscript schedule does not behave like one constant unweighted raw-prefix quantile: its empirical percentile is about 71\\% at depth 1, about 93--95\\% at depths 2 and 3, and above 98\\% at depths 4 through 6. A raw-prefix quantile rule is therefore useful as a transparent calibration device, but by itself it becomes too permissive at deep depths unless the selected quantile is very high or an additional minimum floor/refinement rule is imposed.",
  "",
  "For manuscript framing, the rule can be stated as follows: for a prespecified quantile $q$, define $n_0(d)=\\lceil Q_q\\{m_{d1},\\ldots,m_{dK_d}\\}\\rceil$, where $m_{dk}$ is the number of samples in raw depth-$d$ dominance lineage set $k$ after the fixed AGP discovery filters and before absorb reassignment. The absorb-policy hierarchy is then built using the resulting $n_0(1),\\ldots,n_0(6)$ schedule. The primary analysis should choose one $q$ before outcome testing; sensitivity analysis should repeat the dCST construction and downstream association summaries over a grid such as $q\\in\\{0.75,0.90,0.95,0.975,0.99\\}$, with an explicit statement if a lower bound on $n_0(d)$ is also imposed at deeper levels.",
  "",
  "\\section*{Generated files}",
  "\\begin{itemize}",
  "\\item \\texttt{tables/agp\\_raw\\_lineage\\_size\\_summary.tsv}: depth-specific distribution summary.",
  "\\item \\texttt{tables/agp\\_quantile\\_n0\\_candidates.tsv}: long-format quantile threshold table.",
  "\\item \\texttt{tables/agp\\_raw\\_lineage\\_set\\_sizes\\_by\\_depth.tsv}: all raw lineage-set sizes.",
  "\\item \\texttt{figures/agp\\_lineage\\_size\\_histograms\\_depth1\\_to6.png}: combined histogram figure.",
  "\\item \\texttt{figures/agp\\_lineage\\_size\\_histogram\\_depth*.png}: one histogram per depth.",
  "\\end{itemize}",
  "",
  "\\end{document}"
)
writeLines(report_lines, report_tex)
messagef("Wrote %s", report_tex)
