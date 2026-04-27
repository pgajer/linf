#!/usr/bin/env Rscript
#
# build_calprotectin_context_models.R
#
# Manuscript-ready Halfvarson 2017 Crohn calprotectin follow-up.
#
# This refactored script replaces the earlier prototype context-model workflow
# with a validation-branch-aware pipeline that stays aligned with the sample
# universe used in the current manuscript.
#
# Manuscript-facing analysis:
#   Pragmatic high-inflammatory-burden split at calprotectin >= 250 ug/g
#   plus a threshold sweep across 100, 150, 200, and 250 ug/g
#
# Optional second stage:
#   Mechanistic context-aware modeling can be enabled explicitly via
#   DCST_CALPRO_RUN_CONTEXT_MODEL=1. It is not part of the default manuscript
#   pipeline and is intended for exploratory review before any manuscript
#   integration.

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L || all(is.na(a))) b else a
}

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1] %||% "")
script_path <- if (nzchar(script_arg)) {
  normalizePath(script_arg, winslash = "/", mustWork = TRUE)
} else {
  normalizePath("papers/gut-dcst/scripts/build_calprotectin_context_models.R",
                winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_path)
paper_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
current_projects_root <- normalizePath(file.path(paper_root, "..", "..", ".."),
                                       winslash = "/", mustWork = TRUE)
gut_root <- file.path(current_projects_root, "gut_microbiome")
tables_dir <- file.path(paper_root, "assets", "tables")
figures_dir <- file.path(paper_root, "assets", "figures")
dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(figures_dir, recursive = TRUE, showWarnings = FALSE)

counts_path <- file.path(
  gut_root, "outputs", "validation_ingest", "external_harmonized",
  "halfvarson_2017", "final_matrix", "validation_counts.tsv.gz"
)
validation_meta_path <- file.path(
  gut_root, "outputs", "validation_ingest", "external_harmonized",
  "halfvarson_2017", "validation_metadata_table.tsv"
)
clinical_supp_path <- file.path(
  gut_root, "data", "validation_cohorts", "external", "halfvarson_2017",
  "metadata", "halfvarson_supplementary_dataset1_parsed.tsv"
)
rebuilt_assignments_path <- file.path(
  gut_root, "outputs", "dcst_validation", "absorb", "halfvarson_2017",
  "halfvarson_2017__Crohn_vs_Healthy_dcst_assignments.csv"
)
agp_assignments_path <- file.path(
  gut_root, "outputs", "dcst_validation", "frozen_agp_silva_local_qza_2026-04-26_n0_50_25_25_25", "halfvarson_2017",
  "halfvarson_2017__Crohn_vs_Healthy_dcst_assignments.csv"
)

output_summary_md <- file.path(
  tables_dir, "TABLE_S13H_S13I_halfvarson_calprotectin_threshold_summary.md"
)
output_threshold_tsv <- file.path(
  tables_dir, "TABLE_S13H_halfvarson_calprotectin_threshold250_sensitivity.tsv"
)
output_threshold_tex <- file.path(
  tables_dir, "TABLE_S13H_halfvarson_calprotectin_threshold250_sensitivity.tex"
)
output_threshold_sweep_tsv <- file.path(
  tables_dir, "TABLE_S13I_halfvarson_calprotectin_threshold_sweep.tsv"
)
output_threshold_sweep_tex <- file.path(
  tables_dir, "TABLE_S13I_halfvarson_calprotectin_threshold_sweep.tex"
)
output_threshold_sweep_fig <- file.path(
  figures_dir, "FIGURE_S2_halfvarson_calprotectin_threshold_sweep.png"
)

cd_codes <- c("CCD", "CC", "ICD_r", "ICD_nr", "LC")
depths <- 1:4
threshold_grid <- c(100L, 150L, 200L, 250L)
min_label_count <- suppressWarnings(as.integer(Sys.getenv("DCST_CALPRO_MIN_LABEL_COUNT", "5")))
if (!is.finite(min_label_count) || min_label_count < 2L) min_label_count <- 5L
run_context_model <- tolower(Sys.getenv("DCST_CALPRO_RUN_CONTEXT_MODEL", "0")) %in% c("1", "true", "yes")
context_eps <- suppressWarnings(as.numeric(Sys.getenv("DCST_CALPRO_CONTEXT_EPS", "1e-6")))
if (!is.finite(context_eps) || context_eps <= 0) context_eps <- 1e-6

output_context_feature_map_tsv <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_feature_map.tsv"
)
output_context_feature_map_tex <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_feature_map.tex"
)
output_context_coeff_tsv <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_coefficients.tsv"
)
output_context_coeff_tex <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_coefficients.tex"
)
output_context_slopes_tsv <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_simple_slopes.tsv"
)
output_context_slopes_tex <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_simple_slopes.tex"
)
output_context_corr_tsv <- file.path(
  tables_dir, "TABLE_EXPLORATORY_halfvarson_calprotectin_context_modelC_predictor_correlations.tsv"
)
output_context_partial_fig <- file.path(
  figures_dir, "FIGURE_EXPLORATORY_halfvarson_calprotectin_context_modelC_partial_effects.png"
)
output_context_coef_fig <- file.path(
  figures_dir, "FIGURE_EXPLORATORY_halfvarson_calprotectin_context_modelC_coefficients.png"
)
output_context_summary_md <- file.path(
  tables_dir, "REPORT_EXPLORATORY_halfvarson_calprotectin_context_modelC_summary.md"
)

stop_if_missing <- function(path) {
  if (!file.exists(path)) {
    stop("Required input file is missing: ", path, call. = FALSE)
  }
}

for (path in c(
  counts_path,
  validation_meta_path,
  clinical_supp_path,
  rebuilt_assignments_path,
  agp_assignments_path
)) {
  stop_if_missing(path)
}

if (run_context_model && !requireNamespace("lme4", quietly = TRUE)) {
  stop("Package 'lme4' is required when running the optional context-aware calprotectin model.",
       call. = FALSE)
}

read_delim_auto <- function(path, sep = "\t", select = NULL) {
  if (grepl("\\.gz$", path)) {
    cmd <- sprintf("gzip -dc %s", shQuote(path))
    fread(cmd = cmd, sep = sep, select = select, data.table = FALSE)
  } else {
    fread(path, sep = sep, select = select, data.table = FALSE)
  }
}

read_counts_header <- function(path) {
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "rt") else file(path, open = "rt")
  on.exit(close(con), add = TRUE)
  header <- readLines(con, n = 1L, warn = FALSE)
  if (!length(header) || !nzchar(header)) {
    stop("Unable to read counts header from: ", path, call. = FALSE)
  }
  strsplit(header, "\t", fixed = TRUE)[[1]][-1]
}

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#{}$])", "\\\\\\1", x, perl = TRUE)
  x
}

fmt_num <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf(paste0("%.", digits, "f"), x)
}

format_label_tex <- function(label) {
  if (is.na(label) || !nzchar(label)) return("NA")
  parts <- strsplit(as.character(label), "__", fixed = TRUE)[[1]]
  paste(vapply(parts, tex_escape, character(1L)), collapse = "\\dcstsep ")
}

fmt_int <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

fmt_prob <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  if (x < 1e-4) return(formatC(x, format = "e", digits = 2))
  sprintf("%.4f", x)
}

fmt_or <- function(x) {
  if (is.na(x) || !is.finite(x)) return("NA")
  sprintf("%.3f", x)
}

sanitize_message <- function(x) {
  x <- as.character(x %||% "")
  x <- gsub("[\r\n\t]+", " ", x)
  trimws(x)
}

load_counts_universe <- function() {
  counts <- read_delim_auto(counts_path, select = "run_accession")
  unique(as.character(counts$run_accession))
}

load_halfvarson_clinical <- function() {
  validation_meta <- read_delim_auto(validation_meta_path)
  clinical_meta <- read_delim_auto(clinical_supp_path)

  validation_meta$participant_id_hint <- as.character(validation_meta$participant_id_hint)
  validation_meta$timepoint_hint <- suppressWarnings(as.numeric(validation_meta$timepoint_hint))
  clinical_meta$host_subject_id <- as.character(clinical_meta$host_subject_id)
  clinical_meta$timepoint <- suppressWarnings(as.numeric(clinical_meta$timepoint))

  merged <- merge(
    validation_meta,
    clinical_meta,
    by.x = c("participant_id_hint", "timepoint_hint"),
    by.y = c("host_subject_id", "timepoint"),
    all.x = TRUE
  )
  merged$run_accession <- as.character(merged$run_accession)
  merged$participant_id_hint <- as.character(merged$participant_id_hint)
  merged$calprotectin_num <- suppressWarnings(as.numeric(merged$calprotectin))
  merged
}

load_assignments <- function(mode) {
  mode <- match.arg(mode, c("rebuilt", "agp_transfer"))
  path <- if (mode == "rebuilt") rebuilt_assignments_path else agp_assignments_path
  dt <- read_delim_auto(path, sep = ",")
  needed <- c(
    "run_accession",
    if (mode == "rebuilt") sprintf("depth%d_absorb", depths) else sprintf("depth%d_frozen_agp_absorb", depths)
  )
  missing_cols <- setdiff(needed, names(dt))
  if (length(missing_cols)) {
    stop("Missing expected columns in ", path, ": ", paste(missing_cols, collapse = ", "),
         call. = FALSE)
  }
  dt[, needed, drop = FALSE]
}

assemble_base_frame <- function() {
  counts_runs <- load_counts_universe()
  clinical <- load_halfvarson_clinical()
  clinical <- clinical[
    clinical$run_accession %in% counts_runs &
      clinical$disease_subtype_hint %in% cd_codes &
      is.finite(clinical$calprotectin_num) &
      clinical$calprotectin_num > 0,
    , drop = FALSE
  ]
  clinical$subject_id <- as.character(clinical$participant_id_hint)
  clinical$timepoint_ord <- suppressWarnings(as.numeric(clinical$timepoint_hint))
  clinical$log10_cal <- log10(clinical$calprotectin_num)
  clinical$high_cal_250 <- as.integer(clinical$calprotectin_num >= 250)
  clinical <- clinical[order(clinical$subject_id, clinical$timepoint_ord, clinical$run_accession), , drop = FALSE]
  clinical
}

prepare_mode_frame <- function(base_df, mode) {
  assign <- load_assignments(mode)
  merged <- merge(base_df, assign, by = "run_accession", all.x = FALSE, all.y = FALSE)
  merged <- merged[order(merged$subject_id, merged$timepoint_ord, merged$run_accession), , drop = FALSE]
  merged
}

safe_glmer_label <- function(df, outcome_col) {
  fit_df <- df[, c(outcome_col, "log10_cal", "subject_id"), drop = FALSE]
  names(fit_df)[1] <- "is_label"
  fit_df$is_label <- as.integer(fit_df$is_label)
  fit_df$subject_id <- factor(fit_df$subject_id)
  warnings <- character()
  fit <- withCallingHandlers(
    suppressMessages(
      tryCatch(
        lme4::glmer(
          is_label ~ log10_cal + (1 | subject_id),
          data = fit_df,
          family = binomial(),
          control = lme4::glmerControl(
            optimizer = "bobyqa",
            optCtrl = list(maxfun = 2e5)
          )
        ),
        error = function(e) e
      )
    ),
    warning = function(w) {
      warnings <<- c(warnings, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  list(fit = fit, warnings = unique(warnings))
}

fit_threshold_free_panel <- function(df, label_col, depth, mode_display) {
  panel <- df[!is.na(df[[label_col]]) & nzchar(df[[label_col]]), , drop = FALSE]
  if (!nrow(panel)) return(data.frame())
  labels <- sort(unique(panel[[label_col]]))
  rows <- vector("list", length(labels))
  row_idx <- 0L
  for (label in labels) {
    is_label <- panel[[label_col]] == label
    n_label <- sum(is_label, na.rm = TRUE)
    n_other <- sum(!is_label, na.rm = TRUE)
    if (n_label < min_label_count || n_other < min_label_count) next
    positive_subjects <- length(unique(panel$subject_id[is_label]))
    negative_subjects <- length(unique(panel$subject_id[!is_label]))
    if (positive_subjects < 3L || negative_subjects < 3L) next

    panel$is_label <- is_label
    fit_obj <- safe_glmer_label(panel, "is_label")
    if (inherits(fit_obj$fit, "error")) {
      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        mode = mode_display,
        depth = depth,
        label = label,
        n_samples = nrow(panel),
        n_subjects = length(unique(panel$subject_id)),
        n_label = n_label,
        n_other = n_other,
        beta_log10 = NA_real_,
        se_log10 = NA_real_,
        or_per_doubling = NA_real_,
        ci_low = NA_real_,
        ci_high = NA_real_,
        p_value = NA_real_,
        fit_status = "error",
        fit_message = sanitize_message(conditionMessage(fit_obj$fit)),
        stringsAsFactors = FALSE
      )
      next
    }

    coef_tab <- coef(summary(fit_obj$fit))
    if (!"log10_cal" %in% rownames(coef_tab)) {
      next
    }
    beta <- coef_tab["log10_cal", "Estimate"]
    se <- coef_tab["log10_cal", "Std. Error"]
    p_value <- coef_tab["log10_cal", "Pr(>|z|)"]
    transform_factor <- log10(2)
    or_per_doubling <- exp(beta * transform_factor)
    ci_low <- exp((beta - 1.96 * se) * transform_factor)
    ci_high <- exp((beta + 1.96 * se) * transform_factor)
    warning_text <- sanitize_message(paste(fit_obj$warnings, collapse = " | "))
    status <- if (nzchar(warning_text)) {
      "warning"
    } else if (lme4::isSingular(fit_obj$fit, tol = 1e-4)) {
      "singular"
    } else {
      "ok"
    }
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.frame(
      mode = mode_display,
      depth = depth,
      label = label,
      n_samples = nrow(panel),
      n_subjects = length(unique(panel$subject_id)),
      n_label = n_label,
      n_other = n_other,
      beta_log10 = beta,
      se_log10 = se,
      or_per_doubling = or_per_doubling,
      ci_low = ci_low,
      ci_high = ci_high,
      p_value = p_value,
      fit_status = status,
      fit_message = warning_text,
      stringsAsFactors = FALSE
    )
  }
  if (row_idx == 0L) return(data.frame())
  out <- rbindlist(rows[seq_len(row_idx)], fill = TRUE)
  if (!nrow(out)) return(data.frame())
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$q_value, out$p_value, -out$n_label, out$label), , drop = FALSE]
  rownames(out) <- NULL
  out
}

fit_threshold_panel <- function(df, label_col, depth, mode_display) {
  panel <- df[!is.na(df[[label_col]]) & nzchar(df[[label_col]]), , drop = FALSE]
  if (!nrow(panel)) return(data.frame())
  labels <- sort(unique(panel[[label_col]]))
  rows <- vector("list", length(labels))
  row_idx <- 0L
  for (label in labels) {
    is_label <- panel[[label_col]] == label
    n_label <- sum(is_label, na.rm = TRUE)
    n_other <- sum(!is_label, na.rm = TRUE)
    if (n_label < min_label_count || n_other < min_label_count) next
    a <- sum(panel$high_cal_250 == 1L & is_label, na.rm = TRUE)
    b <- sum(panel$high_cal_250 == 0L & is_label, na.rm = TRUE)
    c <- sum(panel$high_cal_250 == 1L & !is_label, na.rm = TRUE)
    d <- sum(panel$high_cal_250 == 0L & !is_label, na.rm = TRUE)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.frame(
      mode = mode_display,
      depth = depth,
      label = label,
      n_samples = nrow(panel),
      n_subjects = length(unique(panel$subject_id)),
      high_cal_samples = sum(panel$high_cal_250 == 1L, na.rm = TRUE),
      lower_cal_samples = sum(panel$high_cal_250 == 0L, na.rm = TRUE),
      n_high_label = a,
      n_lower_label = b,
      n_high_other = c,
      n_lower_other = d,
      or_value = unname(ft$estimate),
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
  if (row_idx == 0L) return(data.frame())
  out <- rbindlist(rows[seq_len(row_idx)], fill = TRUE)
  if (!nrow(out)) return(data.frame())
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$q_value, out$p_value, -out$n_high_label, out$label), , drop = FALSE]
  rownames(out) <- NULL
  out
}

fit_threshold_panel_at_cutoff <- function(df, label_col, depth, mode_display, threshold) {
  panel <- df[!is.na(df[[label_col]]) & nzchar(df[[label_col]]), , drop = FALSE]
  if (!nrow(panel)) return(data.frame())
  labels <- sort(unique(panel[[label_col]]))
  high_flag <- as.integer(panel$calprotectin_num >= threshold)
  rows <- vector("list", length(labels))
  row_idx <- 0L
  for (label in labels) {
    is_label <- panel[[label_col]] == label
    n_label <- sum(is_label, na.rm = TRUE)
    n_other <- sum(!is_label, na.rm = TRUE)
    if (n_label < min_label_count || n_other < min_label_count) next
    a <- sum(high_flag == 1L & is_label, na.rm = TRUE)
    b <- sum(high_flag == 0L & is_label, na.rm = TRUE)
    c <- sum(high_flag == 1L & !is_label, na.rm = TRUE)
    d <- sum(high_flag == 0L & !is_label, na.rm = TRUE)
    ft <- fisher.test(matrix(c(a, b, c, d), nrow = 2))
    row_idx <- row_idx + 1L
    rows[[row_idx]] <- data.frame(
      threshold = as.integer(threshold),
      mode = mode_display,
      depth = depth,
      label = label,
      n_samples = nrow(panel),
      n_subjects = length(unique(panel$subject_id)),
      high_cal_samples = sum(high_flag == 1L, na.rm = TRUE),
      lower_cal_samples = sum(high_flag == 0L, na.rm = TRUE),
      n_high_label = a,
      n_lower_label = b,
      n_high_other = c,
      n_lower_other = d,
      or_value = unname(ft$estimate),
      ci_low = ft$conf.int[1],
      ci_high = ft$conf.int[2],
      p_value = ft$p.value,
      stringsAsFactors = FALSE
    )
  }
  if (row_idx == 0L) return(data.frame())
  out <- rbindlist(rows[seq_len(row_idx)], fill = TRUE)
  if (!nrow(out)) return(data.frame())
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out <- out[order(out$q_value, out$p_value, -out$n_high_label, out$label), , drop = FALSE]
  rownames(out) <- NULL
  out
}

panel_summary_rows <- function(panel_results, analysis_type) {
  if (!nrow(panel_results)) return(data.frame())
  preferred <- if ("fit_status" %in% names(panel_results)) {
    panel_results[panel_results$fit_status %in% c("ok", "singular"), , drop = FALSE]
  } else {
    panel_results[0, , drop = FALSE]
  }
  top <- if (nrow(preferred)) preferred[1, , drop = FALSE] else panel_results[1, , drop = FALSE]
  top$analysis_type <- analysis_type
  top
}

extract_first_subject_subset <- function(df) {
  idx <- !duplicated(df$subject_id)
  df[idx, , drop = FALSE]
}

build_subject_sensitivity <- function(mode_frames, threshold_free_summary) {
  rows <- list()
  row_idx <- 0L
  for (mode_name in names(mode_frames)) {
    frame <- mode_frames[[mode_name]]
    mode_display <- unique(frame$mode_display)[1]
    for (depth_value in depths) {
      sample_panel <- threshold_free_summary[
        threshold_free_summary$mode == mode_display &
          threshold_free_summary$depth == depth_value,
        , drop = FALSE
      ]
      if (!nrow(sample_panel)) next
      top_sample <- sample_panel[1, , drop = FALSE]
      label_col <- unique(sample_panel$label_col)[1]

      first_df <- extract_first_subject_subset(frame)
      first_results <- fit_threshold_free_panel(first_df, label_col = label_col, depth = depth_value,
                                                mode_display = mode_display)
      first_match <- first_results[first_results$label == top_sample$label, , drop = FALSE]
      row_idx <- row_idx + 1L
      rows[[row_idx]] <- data.frame(
        mode = mode_display,
        depth = depth_value,
        label = top_sample$label,
        sample_or_per_doubling = top_sample$or_per_doubling,
        sample_ci_low = top_sample$ci_low,
        sample_ci_high = top_sample$ci_high,
        sample_q = top_sample$q_value,
        first_subject_samples = nrow(first_df),
        first_subjects = length(unique(first_df$subject_id)),
        first_subject_or_per_doubling = if (nrow(first_match)) first_match$or_per_doubling[[1]] else NA_real_,
        first_subject_ci_low = if (nrow(first_match)) first_match$ci_low[[1]] else NA_real_,
        first_subject_ci_high = if (nrow(first_match)) first_match$ci_high[[1]] else NA_real_,
        first_subject_q = if (nrow(first_match)) first_match$q_value[[1]] else NA_real_,
        stringsAsFactors = FALSE
      )
    }
  }
  if (!row_idx) return(data.frame())
  rbindlist(rows, fill = TRUE)
}

write_threshold_free_tex <- function(summary_df) {
  lines <- c(
    "\\noindent\\textbf{Supplementary Table S13F.} Threshold-free Halfvarson Crohn calprotectin follow-up. Within Crohn-disease samples carrying finite positive calprotectin values, each rebuilt or AGP-transfer dCST label was tested with a mixed-effects logistic model of label membership on continuous log10-calprotectin with subject random intercepts. The table reports the top q-ranked label at each depth and mode. Odds ratios are expressed per doubling of calprotectin rather than per arbitrary threshold split.",
    "",
    "\\begin{center}",
    "\\scriptsize",
    "\\begin{tabularx}{\\textwidth}{@{} l r Y r r Y r r @{}}",
    "\\toprule",
    "Mode & Depth & Top dCST & Samples & Subjects & OR per doubling (95\\% CI) & \\(p\\) & \\(q\\) \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    lines <- c(lines, sprintf(
      "%s & %d & %s & %s & %s & %s [%s, %s] & %s & %s \\\\",
      tex_escape(row$mode),
      as.integer(row$depth),
      format_label_tex(row$label),
      fmt_int(row$n_samples),
      fmt_int(row$n_subjects),
      fmt_or(row$or_per_doubling),
      fmt_or(row$ci_low),
      fmt_or(row$ci_high),
      fmt_prob(row$p_value),
      fmt_prob(row$q_value)
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\end{center}",
    "\\bigskip",
    ""
  )
  writeLines(lines, output_summary_tex)
}

write_subject_sensitivity_tex <- function(subject_df) {
  lines <- c(
    "\\noindent\\textbf{Supplementary Table S13G.} One-sample-per-subject sensitivity for the threshold-free Halfvarson Crohn calprotectin follow-up. For each mode and depth, the sample-level top label from Supplementary Table S13F was re-evaluated after retaining only the first available sample per subject. Odds ratios are expressed per doubling of calprotectin. Cells are reported as NA when the sample-level top label no longer satisfied the minimum-count requirement in the first-sample-per-subject subset.",
    "",
    "\\begin{center}",
    "\\scriptsize",
    "\\begin{tabularx}{\\textwidth}{@{} l r Y Y Y @{}}",
    "\\toprule",
    "Mode & Depth & Top dCST & Sample-level OR [95\\% CI], \\(q\\) & First-sample-per-subject OR [95\\% CI], \\(q\\) \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(subject_df))) {
    row <- subject_df[i, ]
    lines <- c(lines, sprintf(
      "%s & %d & %s & %s [%s, %s], %s & %s [%s, %s], %s \\\\",
      tex_escape(row$mode),
      as.integer(row$depth),
      format_label_tex(row$label),
      fmt_or(row$sample_or_per_doubling),
      fmt_or(row$sample_ci_low),
      fmt_or(row$sample_ci_high),
      fmt_prob(row$sample_q),
      fmt_or(row$first_subject_or_per_doubling),
      fmt_or(row$first_subject_ci_low),
      fmt_or(row$first_subject_ci_high),
      fmt_prob(row$first_subject_q)
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\end{center}",
    "\\bigskip",
    ""
  )
  writeLines(lines, output_subject_tex)
}

write_threshold_tex <- function(summary_df) {
  lines <- c(
    "\\noindent\\textbf{Supplementary Table S13H.} Thresholded sensitivity for the Halfvarson Crohn calprotectin follow-up. Using the joined Crohn-disease sample universe from the Halfvarson validation matrix, calprotectin was dichotomized at 250~\\textmu g/g as a pragmatic high-inflammatory-burden contrast and each label panel was re-evaluated by Fisher's exact test. The table reports the top q-ranked label at each depth and mode.",
    "",
    "\\begin{center}",
    "\\scriptsize",
    "\\begin{tabularx}{\\textwidth}{@{} l r Y r r Y r r @{}}",
    "\\toprule",
    "Mode & Depth & Top dCST & High/low samples & Subjects & OR (95\\% CI) & \\(p\\) & \\(q\\) \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    lines <- c(lines, sprintf(
      "%s & %d & %s & %s/%s & %s & %s [%s, %s] & %s & %s \\\\",
      tex_escape(row$mode),
      as.integer(row$depth),
      format_label_tex(row$label),
      fmt_int(row$high_cal_samples),
      fmt_int(row$lower_cal_samples),
      fmt_int(row$n_subjects),
      fmt_or(row$or_value),
      fmt_or(row$ci_low),
      fmt_or(row$ci_high),
      fmt_prob(row$p_value),
      fmt_prob(row$q_value)
    ))
  }
  lines <- c(
    lines,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\end{center}",
    "\\bigskip",
    ""
  )
  writeLines(lines, output_threshold_tex)
}

write_threshold_sweep_tex <- function(summary_df) {
  lines <- c(
    "\\noindent\\textbf{Supplementary Table S13I.} Threshold-sweep sensitivity for the Halfvarson Crohn calprotectin follow-up. Using the same Crohn-disease sample universe as Supplementary Table S13H, calprotectin was dichotomized at 100, 150, 200, and 250~\\textmu g/g and the top q-ranked label at each depth and mode was recorded. This sweep is intended to show whether the thresholded signal is concentrated only at clearly elevated calprotectin values or remains stable across lower clinical cut-offs.",
    "",
    "\\scriptsize",
    "\\setlength{\\LTleft}{0pt}",
    "\\setlength{\\LTright}{0pt}",
    "\\begin{longtable}{@{} r l r p{0.30\\textwidth} r r p{0.20\\textwidth} r r @{}}",
    "\\toprule",
    "Threshold & Mode & Depth & Top dCST & High/low samples & Subjects & OR (95\\% CI) & \\(p\\) & \\(q\\) \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    "Threshold & Mode & Depth & Top dCST & High/low samples & Subjects & OR (95\\% CI) & \\(p\\) & \\(q\\) \\\\",
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{9}{r}{\\emph{Continued on next page}} \\\\",
    "\\midrule",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot"
  )
  for (i in seq_len(nrow(summary_df))) {
    row <- summary_df[i, ]
    lines <- c(lines, sprintf(
      "%d & %s & %d & %s & %s/%s & %s & %s [%s, %s] & %s & %s \\\\",
      as.integer(row$threshold),
      tex_escape(row$mode),
      as.integer(row$depth),
      format_label_tex(row$label),
      fmt_int(row$high_cal_samples),
      fmt_int(row$lower_cal_samples),
      fmt_int(row$n_subjects),
      fmt_or(row$or_value),
      fmt_or(row$ci_low),
      fmt_or(row$ci_high),
      fmt_prob(row$p_value),
      fmt_prob(row$q_value)
    ))
  }
  lines <- c(lines, "\\end{longtable}", "\\normalsize")
  writeLines(lines, output_threshold_sweep_tex)
}

plot_threshold_sweep_figure <- function(summary_df) {
  if (!nrow(summary_df)) {
    return(invisible(NULL))
  }

  mode_order <- c("Rebuilt", "AGP transfer")
  depth_colors <- c(
    "1" = "#4c78a8",
    "2" = "#f58518",
    "3" = "#54a24b",
    "4" = "#b279a2"
  )
  depth_pch <- c("1" = 16, "2" = 17, "3" = 15, "4" = 18)
  threshold_values <- sort(unique(summary_df$threshold))
  q_line <- -log10(0.05)
  y_max <- max(-log10(summary_df$q_value), q_line, na.rm = TRUE)
  y_max <- max(y_max + 0.15, 1.6)

  grDevices::png(
    filename = output_threshold_sweep_fig,
    width = 2200,
    height = 1100,
    res = 220
  )
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)

  par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3.4, 1.2), oma = c(0, 0, 1.3, 0))

  for (mode_name in mode_order) {
    panel_df <- summary_df[summary_df$mode == mode_name, , drop = FALSE]
    plot(
      NA,
      xlim = range(threshold_values),
      ylim = c(0, y_max),
      xaxt = "n",
      xlab = "Calprotectin threshold (µg/g)",
      ylab = expression(-log[10](q)),
      main = mode_name,
      bty = "l"
    )
    axis(1, at = threshold_values, labels = threshold_values)
    abline(h = q_line, lty = 2, lwd = 1.3, col = "#7a7a7a")
    grid(nx = NA, ny = NULL, col = "#ececec", lty = 1)

    for (depth_value in depths) {
      depth_key <- as.character(depth_value)
      depth_df <- panel_df[panel_df$depth == depth_value, , drop = FALSE]
      if (!nrow(depth_df)) next
      depth_df <- depth_df[order(depth_df$threshold), , drop = FALSE]
      y_vals <- -log10(depth_df$q_value)
      lines(
        depth_df$threshold,
        y_vals,
        col = depth_colors[[depth_key]],
        lwd = 2.2,
        type = "b",
        pch = depth_pch[[depth_key]],
        cex = 1.0
      )
    }

    if (mode_name == "Rebuilt") {
      legend(
        "topright",
        legend = paste("Depth", depths),
        col = unname(depth_colors[as.character(depths)]),
        pch = unname(depth_pch[as.character(depths)]),
        lwd = 2.2,
        pt.cex = 1,
        bty = "n",
        inset = 0.01
      )
    }
  }

  mtext(
    "Halfvarson Crohn calprotectin threshold sweep",
    outer = TRUE,
    cex = 1.15,
    font = 2
  )
}

match_context_features <- function(header) {
  pick <- function(patterns) {
    unique(unlist(lapply(
      patterns,
      function(pattern) header[grepl(pattern, header, perl = TRUE, ignore.case = TRUE)]
    )))
  }

  fp_feature <- header[grepl("^Faecalibacterium prausnitzii$", header, perl = TRUE, ignore.case = TRUE)]
  if (!length(fp_feature)) {
    fp_feature <- header[grepl("^Faecalibacterium\\b", header, perl = TRUE, ignore.case = TRUE)][1]
  }
  if (!length(fp_feature) || is.na(fp_feature) || !nzchar(fp_feature)) {
    stop("Could not identify a Faecalibacterium prausnitzii feature in the Halfvarson counts matrix.",
         call. = FALSE)
  }

  feature_map <- list(
    Fp_feature = fp_feature[[1]],
    cross_feeding = list(
      requested = c("Bifidobacterium", "Blautia", "Roseburia"),
      matched = pick(c("^Bifidobacterium\\b", "^Blautia\\b", "^Roseburia\\b"))
    ),
    oxy_tolerators = list(
      requested = c("Escherichia-Shigella", "Enterobacteriaceae", "Klebsiella"),
      matched = pick(c("^Escherichia-Shigella\\b", "^Enterobacteriaceae\\b", "^Klebsiella\\b"))
    ),
    bacteroides_bg = list(
      requested = c("Bacteroides"),
      matched = pick(c("^Bacteroides\\b", "^\\[Bacteroides\\]"))
    ),
    mucin_context = list(
      requested = c("Akkermansia", "Ruminococcus_torques"),
      matched = pick(c("^Akkermansia\\b", "^Ruminococcus[ _]torques\\b"))
    )
  )
  feature_map
}

build_context_model_frame <- function(base_df) {
  header <- read_counts_header(counts_path)
  feature_map <- match_context_features(header)
  selected_features <- unique(c(
    feature_map$Fp_feature,
    unlist(lapply(feature_map[c("cross_feeding", "oxy_tolerators", "bacteroides_bg", "mucin_context")],
                  `[[`, "matched"), use.names = FALSE)
  ))
  counts_subset <- read_delim_auto(
    counts_path,
    select = c("run_accession", selected_features)
  )
  merged <- merge(base_df, counts_subset, by = "run_accession", all.x = FALSE, all.y = FALSE)
  merged <- merged[order(merged$subject_id, merged$timepoint_ord, merged$run_accession), , drop = FALSE]
  duplicate_runs_removed <- sum(duplicated(merged$run_accession))
  if (duplicate_runs_removed > 0L) {
    merged <- merged[!duplicated(merged$run_accession), , drop = FALSE]
  }

  numeric_mat <- as.matrix(merged[, selected_features, drop = FALSE])
  storage.mode(numeric_mat) <- "numeric"
  lib_size <- rowSums(numeric_mat, na.rm = TRUE)
  lib_size[!is.finite(lib_size) | lib_size <= 0] <- NA_real_
  rel_mat <- numeric_mat / lib_size
  colnames(rel_mat) <- selected_features

  score_sum <- function(features) {
    if (!length(features)) {
      rep(0, nrow(rel_mat))
    } else {
      rowSums(rel_mat[, features, drop = FALSE], na.rm = TRUE)
    }
  }

  merged$Fp_rel <- rel_mat[, feature_map$Fp_feature]
  merged$cross_feeding_rel <- score_sum(feature_map$cross_feeding$matched)
  merged$oxy_tolerators_rel <- score_sum(feature_map$oxy_tolerators$matched)
  merged$bacteroides_bg_rel <- score_sum(feature_map$bacteroides_bg$matched)
  merged$mucin_context_rel <- score_sum(feature_map$mucin_context$matched)

  merged$log10_Fp <- log10(merged$Fp_rel + context_eps)
  merged$cross_feeding <- log10(merged$cross_feeding_rel + context_eps)
  merged$oxy_tolerators <- log10(merged$oxy_tolerators_rel + context_eps)
  merged$bacteroides_bg <- log10(merged$bacteroides_bg_rel + context_eps)
  merged$mucin_context <- log10(merged$mucin_context_rel + context_eps)

  keep <- complete.cases(merged[, c(
    "log10_cal", "log10_Fp", "cross_feeding", "oxy_tolerators",
    "bacteroides_bg", "mucin_context", "subject_id"
  )])
  merged <- merged[keep, , drop = FALSE]

  attr(merged, "feature_map") <- feature_map
  attr(merged, "duplicate_runs_removed") <- duplicate_runs_removed
  merged
}

fit_context_model <- function(context_df) {
  fit <- lme4::lmer(
    log10_cal ~ log10_Fp * cross_feeding +
      log10_Fp * oxy_tolerators +
      log10_Fp * bacteroides_bg +
      log10_Fp * mucin_context +
      (1 | subject_id),
    data = context_df,
    REML = FALSE
  )
  fit
}

coefficient_p_table <- function(fit) {
  coef_tab <- coef(summary(fit))
  out <- data.frame(
    term = rownames(coef_tab),
    estimate = coef_tab[, "Estimate"],
    std_error = coef_tab[, "Std. Error"],
    statistic = coef_tab[, ncol(coef_tab)],
    stringsAsFactors = FALSE,
    row.names = NULL
  )
  vc <- as.matrix(stats::vcov(fit))
  out$ci_low <- out$estimate - 1.96 * out$std_error
  out$ci_high <- out$estimate + 1.96 * out$std_error
  p_col <- grep("^Pr\\(", colnames(coef_tab), value = TRUE)
  if (length(p_col)) {
    out$p_value <- coef_tab[, p_col[1]]
  } else {
    out$p_value <- 2 * stats::pnorm(-abs(out$statistic))
  }
  out$term_group <- ifelse(grepl(":", out$term, fixed = TRUE), "interaction",
                           ifelse(out$term == "(Intercept)", "intercept", "main_effect"))
  out
}

context_feature_map_table <- function(feature_map) {
  rows <- list(
    data.frame(
      quantity = "Faecalibacterium prausnitzii anchor",
      requested_components = "Faecalibacterium prausnitzii",
      matched_features = feature_map$Fp_feature,
      n_matched = 1L,
      stringsAsFactors = FALSE
    )
  )
  idx <- 1L
  for (nm in c("cross_feeding", "oxy_tolerators", "bacteroides_bg", "mucin_context")) {
    idx <- idx + 1L
    matched <- feature_map[[nm]]$matched
    rows[[idx]] <- data.frame(
      quantity = nm,
      requested_components = paste(feature_map[[nm]]$requested, collapse = ", "),
      matched_features = if (length(matched)) paste(matched, collapse = "; ") else "none in matrix",
      n_matched = length(matched),
      stringsAsFactors = FALSE
    )
  }
  rbindlist(rows, fill = TRUE)
}

context_predictor_correlations <- function(context_df) {
  vars <- c("log10_Fp", "cross_feeding", "oxy_tolerators", "bacteroides_bg", "mucin_context")
  cor_mat <- stats::cor(context_df[, vars, drop = FALSE], use = "pairwise.complete.obs")
  as.data.frame(as.table(cor_mat), stringsAsFactors = FALSE)
}

interaction_term_name <- function(context_var) {
  candidates <- c(
    sprintf("log10_Fp:%s", context_var),
    sprintf("%s:log10_Fp", context_var)
  )
  candidates
}

extract_fixef_value <- function(fixef_vec, term, default = 0) {
  if (term %in% names(fixef_vec)) fixef_vec[[term]] else default
}

simple_slope_rows <- function(fit, context_df) {
  fixef_vec <- lme4::fixef(fit)
  vc <- as.matrix(stats::vcov(fit))
  context_vars <- c("cross_feeding", "oxy_tolerators", "bacteroides_bg", "mucin_context")
  medians <- vapply(context_vars, function(v) stats::median(context_df[[v]], na.rm = TRUE), numeric(1L))
  levels_map <- c("low (Q1)" = 0.25, "median" = 0.5, "high (Q3)" = 0.75)
  rows <- list()
  idx <- 0L

  for (context_var in context_vars) {
    q_vals <- stats::quantile(context_df[[context_var]], probs = unname(levels_map), na.rm = TRUE)
    names(q_vals) <- names(levels_map)
    for (lvl in names(levels_map)) {
      context_values <- medians
      context_values[[context_var]] <- unname(q_vals[[lvl]])
      contrast <- stats::setNames(rep(0, length(fixef_vec)), names(fixef_vec))
      contrast["log10_Fp"] <- 1
      for (ctx_name in context_vars) {
        int_term <- interaction_term_name(ctx_name)
        int_term <- int_term[int_term %in% names(contrast)][1] %||% NA_character_
        if (!is.na(int_term)) {
          contrast[int_term] <- context_values[[ctx_name]]
        }
      }
      slope <- sum(contrast * fixef_vec[names(contrast)])
      slope_se <- sqrt(drop(t(contrast) %*% vc[names(contrast), names(contrast), drop = FALSE] %*% contrast))
      idx <- idx + 1L
      rows[[idx]] <- data.frame(
        context = context_var,
        level = lvl,
        context_value = unname(q_vals[[lvl]]),
        slope_log10Fp = slope,
        slope_ci_low = slope - 1.96 * slope_se,
        slope_ci_high = slope + 1.96 * slope_se,
        stringsAsFactors = FALSE
      )
    }
  }
  rbindlist(rows, fill = TRUE)
}

predict_fixed_context <- function(fixef_vec, log10_Fp, context_values) {
  intercept <- extract_fixef_value(fixef_vec, "(Intercept)")
  mu <- intercept +
    extract_fixef_value(fixef_vec, "log10_Fp") * log10_Fp +
    extract_fixef_value(fixef_vec, "cross_feeding") * context_values[["cross_feeding"]] +
    extract_fixef_value(fixef_vec, "oxy_tolerators") * context_values[["oxy_tolerators"]] +
    extract_fixef_value(fixef_vec, "bacteroides_bg") * context_values[["bacteroides_bg"]] +
    extract_fixef_value(fixef_vec, "mucin_context") * context_values[["mucin_context"]]

  for (ctx_name in names(context_values)) {
    int_term <- interaction_term_name(ctx_name)
    int_term <- int_term[int_term %in% names(fixef_vec)][1] %||% NA_character_
    if (!is.na(int_term)) {
      mu <- mu + extract_fixef_value(fixef_vec, int_term) * log10_Fp * context_values[[ctx_name]]
    }
  }
  mu
}

plot_context_partial_effects <- function(fit, context_df) {
  fixef_vec <- lme4::fixef(fit)
  context_vars <- c("cross_feeding", "oxy_tolerators", "bacteroides_bg", "mucin_context")
  labels <- c(
    cross_feeding = "Cross-feeding score",
    oxy_tolerators = "Oxygen-tolerator score",
    bacteroides_bg = "Bacteroides background",
    mucin_context = "Mucin context"
  )
  medians <- vapply(context_vars, function(v) stats::median(context_df[[v]], na.rm = TRUE), numeric(1L))
  x_grid <- seq(
    stats::quantile(context_df$log10_Fp, 0.05, na.rm = TRUE),
    stats::quantile(context_df$log10_Fp, 0.95, na.rm = TRUE),
    length.out = 100
  )
  level_probs <- c("Q1" = 0.25, "Median" = 0.5, "Q3" = 0.75)
  level_cols <- c("Q1" = "#4c78a8", "Median" = "#f58518", "Q3" = "#54a24b")

  grDevices::png(output_context_partial_fig, width = 1800, height = 1600, res = 220)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  par(mfrow = c(2, 2), mar = c(4.3, 4.6, 3.0, 1.3), oma = c(0, 0, 1.5, 0))

  for (context_var in context_vars) {
    plot(
      context_df$log10_Fp, context_df$log10_cal,
      col = grDevices::rgb(0, 0, 0, 0.12),
      pch = 16, cex = 0.7,
      xlab = expression(log[10](F.~prausnitzii + epsilon)),
      ylab = expression(log[10](calprotectin)),
      main = labels[[context_var]]
    )
    q_vals <- stats::quantile(context_df[[context_var]], probs = unname(level_probs), na.rm = TRUE)
    names(q_vals) <- names(level_probs)
    for (lvl in names(level_probs)) {
      context_values <- medians
      context_values[[context_var]] <- unname(q_vals[[lvl]])
      preds <- vapply(x_grid, function(xv) predict_fixed_context(fixef_vec, xv, context_values), numeric(1L))
      lines(x_grid, preds, col = level_cols[[lvl]], lwd = 2.4)
    }
    if (context_var == "cross_feeding") {
      legend(
        "topleft",
        legend = sprintf("%s (%s)", names(level_probs), vapply(q_vals, fmt_num, character(1L), digits = 2)),
        col = level_cols[names(level_probs)],
        lwd = 2.4,
        bty = "n",
        cex = 0.9
      )
    }
  }
  mtext("Halfvarson Crohn calprotectin context model: partial effects", outer = TRUE, cex = 1.15, font = 2)
}

plot_context_coefficients <- function(coeff_df) {
  keep <- coeff_df[coeff_df$term != "(Intercept)", , drop = FALSE]
  if (!nrow(keep)) return(invisible(NULL))
  keep$display_term <- keep$term
  keep$display_term <- gsub("log10_Fp", "log10(Fp)", keep$display_term, fixed = TRUE)
  keep$display_term <- gsub("cross_feeding", "cross-feeding", keep$display_term, fixed = TRUE)
  keep$display_term <- gsub("oxy_tolerators", "oxygen-tolerators", keep$display_term, fixed = TRUE)
  keep$display_term <- gsub("bacteroides_bg", "Bacteroides background", keep$display_term, fixed = TRUE)
  keep$display_term <- gsub("mucin_context", "mucin context", keep$display_term, fixed = TRUE)
  keep <- keep[order(keep$estimate), , drop = FALSE]

  grDevices::png(output_context_coef_fig, width = 2200, height = 1400, res = 220)
  old_par <- par(no.readonly = TRUE)
  on.exit({
    par(old_par)
    grDevices::dev.off()
  }, add = TRUE)
  par(mar = c(5.0, 18.0, 3.0, 2.0))
  y_pos <- seq_len(nrow(keep))
  plot(
    keep$estimate, y_pos,
    xlim = range(c(keep$ci_low, keep$ci_high), na.rm = TRUE),
    ylim = c(0.5, nrow(keep) + 0.5),
    yaxt = "n",
    xlab = "Coefficient estimate on log10(calprotectin)",
    ylab = "",
    pch = 16,
    col = ifelse(keep$term_group == "interaction", "#d62728", "#4c78a8")
  )
  axis(2, at = y_pos, labels = keep$display_term, las = 1)
  abline(v = 0, lty = 2, col = "#7a7a7a")
  segments(keep$ci_low, y_pos, keep$ci_high, y_pos,
           col = ifelse(keep$term_group == "interaction", "#d62728", "#4c78a8"),
           lwd = 2.0)
  title("Halfvarson Crohn calprotectin context model coefficients", cex.main = 1.3)
}

write_context_feature_map_tex <- function(feature_df) {
  lines <- c(
    "\\noindent\\textbf{Exploratory context-score feature map.} Requested taxa for the mechanistic context model and the actual Halfvarson validation-matrix features matched locally.",
    "",
    "\\begin{center}",
    "\\scriptsize",
    "\\begin{tabularx}{\\textwidth}{@{} l Y Y r @{}}",
    "\\toprule",
    "Quantity & Requested components & Matched matrix features & Count \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(feature_df))) {
    row <- feature_df[i, ]
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s \\\\",
      tex_escape(row$quantity),
      tex_escape(row$requested_components),
      tex_escape(row$matched_features),
      fmt_int(row$n_matched)
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabularx}", "\\end{center}", "")
  writeLines(lines, output_context_feature_map_tex)
}

write_context_coeff_tex <- function(coeff_df) {
  lines <- c(
    "\\noindent\\textbf{Exploratory multivariable context-model coefficients.} Linear mixed-effects model for continuous log10-calprotectin in Halfvarson Crohn samples with subject random intercepts.",
    "",
    "\\begin{center}",
    "\\scriptsize",
    "\\begin{tabularx}{\\textwidth}{@{} l r r r r r @{}}",
    "\\toprule",
    "Term & Estimate & SE & 95\\% CI low & 95\\% CI high & Approx. \\(p\\) \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(coeff_df))) {
    row <- coeff_df[i, ]
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & %s & %s \\\\",
      tex_escape(row$term),
      fmt_num(row$estimate),
      fmt_num(row$std_error),
      fmt_num(row$ci_low),
      fmt_num(row$ci_high),
      fmt_prob(row$p_value)
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabularx}", "\\end{center}", "")
  writeLines(lines, output_context_coeff_tex)
}

write_context_slopes_tex <- function(slopes_df) {
  lines <- c(
    "\\noindent\\textbf{Exploratory simple slopes for Faecalibacterium prausnitzii.} Estimated slope of log10-calprotectin on log10(F. prausnitzii abundance) while varying one context score at a time and holding the other context scores at their sample medians.",
    "",
    "\\begin{center}",
    "\\scriptsize",
    "\\begin{tabularx}{\\textwidth}{@{} l l r r r @{}}",
    "\\toprule",
    "Context & Level & Context value & Slope & 95\\% CI \\\\",
    "\\midrule"
  )
  for (i in seq_len(nrow(slopes_df))) {
    row <- slopes_df[i, ]
    lines <- c(lines, sprintf(
      "%s & %s & %s & %s & [%s, %s] \\\\",
      tex_escape(row$context),
      tex_escape(row$level),
      fmt_num(row$context_value),
      fmt_num(row$slope_log10Fp),
      fmt_num(row$slope_ci_low),
      fmt_num(row$slope_ci_high)
    ))
  }
  lines <- c(lines, "\\bottomrule", "\\end{tabularx}", "\\end{center}", "")
  writeLines(lines, output_context_slopes_tex)
}

write_context_summary_markdown <- function(context_df, feature_df, coeff_df, slopes_df) {
  context_terms <- coeff_df[coeff_df$term_group == "interaction", , drop = FALSE]
  duplicate_runs_removed <- attr(context_df, "duplicate_runs_removed") %||% 0L
  lines <- c(
    "# Exploratory Halfvarson Crohn calprotectin context model",
    "",
    sprintf("- Sample universe: %s Crohn runs across %s subjects", fmt_int(nrow(context_df)),
            fmt_int(length(unique(context_df$subject_id)))),
    sprintf("- Exact duplicated run accessions removed before fitting: %s", fmt_int(duplicate_runs_removed)),
    sprintf("- Epsilon added before log10 transforms: %s", format(context_eps, scientific = TRUE)),
    "",
    "## Actual feature map",
    ""
  )
  for (i in seq_len(nrow(feature_df))) {
    row <- feature_df[i, ]
    lines <- c(lines, sprintf(
      "- %s: matched %s feature(s): %s",
      row$quantity, fmt_int(row$n_matched), row$matched_features
    ))
  }
  lines <- c(lines, "", "## Interaction-term summary", "")
  if (nrow(context_terms)) {
    for (i in seq_len(nrow(context_terms))) {
      row <- context_terms[i, ]
      lines <- c(lines, sprintf(
        "- %s: estimate %s [%s, %s], approx. p = %s",
        row$term, fmt_num(row$estimate), fmt_num(row$ci_low), fmt_num(row$ci_high), fmt_prob(row$p_value)
      ))
    }
  }
  lines <- c(lines, "", "## Simple slopes", "")
  for (i in seq_len(nrow(slopes_df))) {
    row <- slopes_df[i, ]
    lines <- c(lines, sprintf(
      "- %s | %s: slope %s [%s, %s]",
      row$context, row$level, fmt_num(row$slope_log10Fp), fmt_num(row$slope_ci_low), fmt_num(row$slope_ci_high)
    ))
  }
  writeLines(lines, output_context_summary_md)
}

run_context_model_stage <- function(base_df) {
  message("=== Exploratory Halfvarson calprotectin context model (Option C) ===")
  context_df <- build_context_model_frame(base_df)
  feature_df <- context_feature_map_table(attr(context_df, "feature_map"))
  fit <- fit_context_model(context_df)
  coeff_df <- coefficient_p_table(fit)
  slopes_df <- simple_slope_rows(fit, context_df)
  corr_df <- context_predictor_correlations(context_df)

  fwrite(feature_df, output_context_feature_map_tsv, sep = "\t", quote = FALSE, na = "NA")
  fwrite(coeff_df, output_context_coeff_tsv, sep = "\t", quote = FALSE, na = "NA")
  fwrite(slopes_df, output_context_slopes_tsv, sep = "\t", quote = FALSE, na = "NA")
  fwrite(corr_df, output_context_corr_tsv, sep = "\t", quote = FALSE, na = "NA")
  write_context_feature_map_tex(feature_df)
  write_context_coeff_tex(coeff_df)
  write_context_slopes_tex(slopes_df)
  plot_context_partial_effects(fit, context_df)
  plot_context_coefficients(coeff_df)
  write_context_summary_markdown(context_df, feature_df, coeff_df, slopes_df)

  message("Exploratory context-model outputs:")
  message("  ", output_context_feature_map_tsv)
  message("  ", output_context_coeff_tsv)
  message("  ", output_context_slopes_tsv)
  message("  ", output_context_corr_tsv)
  message("  ", output_context_partial_fig)
  message("  ", output_context_coef_fig)
  message("  ", output_context_summary_md)

  list(
    context_df = context_df,
    feature_df = feature_df,
    coeff_df = coeff_df,
    slopes_df = slopes_df,
    corr_df = corr_df,
    fit = fit
  )
}

write_summary_markdown <- function(base_df, threshold_summary, threshold_sweep_summary) {
  repeated_subjects <- sum(table(base_df$subject_id) > 1L)
  lines <- c(
    "# Halfvarson calprotectin threshold sensitivity summary",
    "",
    "## Sample universe",
    sprintf("- Crohn runs with finite positive calprotectin in the Halfvarson validation matrix: %s", fmt_int(nrow(base_df))),
    sprintf("- Unique subjects: %s", fmt_int(length(unique(base_df$subject_id)))),
    sprintf("- Repeated-sample subjects: %s", fmt_int(repeated_subjects)),
    sprintf("- Calprotectin threshold sensitivity split at 250 ug/g: %s high / %s lower visits",
            fmt_int(sum(base_df$high_cal_250 == 1L)), fmt_int(sum(base_df$high_cal_250 == 0L))),
    "",
    "## Outputs",
    sprintf("- %s", basename(output_threshold_tsv)),
    sprintf("- %s", basename(output_threshold_tex)),
    sprintf("- %s", basename(output_threshold_sweep_tsv)),
    sprintf("- %s", basename(output_threshold_sweep_tex)),
    sprintf("- %s", basename(output_threshold_sweep_fig)),
    "",
    "## Notes",
    "- The manuscript-facing calprotectin follow-up is thresholded and sensitivity-oriented.",
    "- The 250 ug/g table is retained as a clinically legible high-inflammatory-burden layer.",
    "- Mechanistic context-aware F. prausnitzii modeling is intentionally left as a second-stage extension rather than part of this primary pipeline.",
    ""
  )
  if (nrow(threshold_summary)) {
    lines <- c(lines, "## Threshold-250 top signals", "")
    for (i in seq_len(nrow(threshold_summary))) {
      row <- threshold_summary[i, ]
      lines <- c(lines, sprintf(
        "- %s depth %d: %s | OR %s [%s, %s], q = %s",
        row$mode,
        as.integer(row$depth),
        row$label,
        fmt_or(row$or_value),
        fmt_or(row$ci_low),
        fmt_or(row$ci_high),
        fmt_prob(row$q_value)
      ))
    }
  }
  if (nrow(threshold_sweep_summary)) {
    lines <- c(lines, "", "## Threshold-sweep top signals", "")
    for (i in seq_len(nrow(threshold_sweep_summary))) {
      row <- threshold_sweep_summary[i, ]
      lines <- c(lines, sprintf(
        "- %d ug/g | %s depth %d: %s | OR %s [%s, %s], q = %s",
        as.integer(row$threshold),
        row$mode,
        as.integer(row$depth),
        row$label,
        fmt_or(row$or_value),
        fmt_or(row$ci_low),
        fmt_or(row$ci_high),
        fmt_prob(row$q_value)
      ))
    }
  }
  writeLines(lines, output_summary_md)
}

main <- function() {
  message("=== Halfvarson calprotectin threshold sensitivity build ===")
  base_df <- assemble_base_frame()
  message("Base frame: ", nrow(base_df), " Crohn samples across ",
          length(unique(base_df$subject_id)), " subjects.")

  mode_frames <- list(
    rebuilt = prepare_mode_frame(base_df, "rebuilt"),
    agp_transfer = prepare_mode_frame(base_df, "agp_transfer")
  )
  mode_frames$rebuilt$mode_display <- "Rebuilt"
  mode_frames$agp_transfer$mode_display <- "AGP transfer"

  threshold_rows <- list()
  threshold_idx <- 0L
  threshold_sweep_rows <- list()
  threshold_sweep_idx <- 0L

  for (mode_name in names(mode_frames)) {
    frame <- mode_frames[[mode_name]]
    mode_display <- unique(frame$mode_display)[1]
    for (depth in depths) {
      label_col <- if (mode_name == "rebuilt") {
        sprintf("depth%d_absorb", depth)
      } else {
        sprintf("depth%d_frozen_agp_absorb", depth)
      }
      th_panel <- fit_threshold_panel(frame, label_col, depth, mode_display)
      if (nrow(th_panel)) {
        th_panel$mode_key <- mode_name
        th_panel$label_col <- label_col
        threshold_idx <- threshold_idx + 1L
        threshold_rows[[threshold_idx]] <- th_panel
      }
      for (threshold_value in threshold_grid) {
        sweep_panel <- fit_threshold_panel_at_cutoff(frame, label_col, depth, mode_display, threshold_value)
        if (nrow(sweep_panel)) {
          sweep_panel$mode_key <- mode_name
          sweep_panel$label_col <- label_col
          threshold_sweep_idx <- threshold_sweep_idx + 1L
          threshold_sweep_rows[[threshold_sweep_idx]] <- sweep_panel
        }
      }
    }
  }

  threshold_master <- if (length(threshold_rows)) {
    rbindlist(threshold_rows, fill = TRUE)
  } else {
    data.frame()
  }
  threshold_sweep_master <- if (length(threshold_sweep_rows)) {
    rbindlist(threshold_sweep_rows, fill = TRUE)
  } else {
    data.frame()
  }

  fwrite(threshold_master, output_threshold_tsv, sep = "\t", quote = FALSE, na = "NA")
  fwrite(threshold_sweep_master, output_threshold_sweep_tsv, sep = "\t", quote = FALSE, na = "NA")

  threshold_summary <- if (nrow(threshold_master)) {
    rbindlist(
      lapply(
        split(threshold_master,
              list(threshold_master$mode, threshold_master$depth), drop = TRUE),
        panel_summary_rows,
        analysis_type = "threshold_250"
      ),
      fill = TRUE
    )
  } else {
    data.frame()
  }
  if (nrow(threshold_summary)) {
    threshold_summary <- threshold_summary[
      order(match(threshold_summary$mode, c("Rebuilt", "AGP transfer")),
            threshold_summary$depth),
      , drop = FALSE
    ]
  }

  threshold_sweep_summary <- if (nrow(threshold_sweep_master)) {
    rbindlist(
      lapply(
        split(
          threshold_sweep_master,
          list(threshold_sweep_master$threshold,
               threshold_sweep_master$mode,
               threshold_sweep_master$depth),
          drop = TRUE
        ),
        panel_summary_rows,
        analysis_type = "threshold_sweep"
      ),
      fill = TRUE
    )
  } else {
    data.frame()
  }
  if (nrow(threshold_sweep_summary)) {
    threshold_sweep_summary <- threshold_sweep_summary[
      order(threshold_sweep_summary$threshold,
            match(threshold_sweep_summary$mode, c("Rebuilt", "AGP transfer")),
            threshold_sweep_summary$depth),
      , drop = FALSE
    ]
  }

  write_threshold_tex(threshold_summary)
  write_threshold_sweep_tex(threshold_sweep_summary)
  plot_threshold_sweep_figure(threshold_sweep_summary)
  write_summary_markdown(base_df, threshold_summary, threshold_sweep_summary)

  message("Wrote:")
  message("  ", output_threshold_tsv)
  message("  ", output_threshold_tex)
  message("  ", output_threshold_sweep_tsv)
  message("  ", output_threshold_sweep_tex)
  message("  ", output_threshold_sweep_fig)
  message("  ", output_summary_md)

  if (run_context_model) {
    run_context_model_stage(base_df)
  }
}

if (sys.nframe() == 0L) {
  invisible(main())
}
