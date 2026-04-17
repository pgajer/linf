#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

`%||%` <- function(a, b) {
  if (is.null(a) || length(a) == 0L || is.na(a)) b else a
}

script_args <- commandArgs(trailingOnly = FALSE)
script_arg <- sub("^--file=", "", grep("^--file=", script_args, value = TRUE)[1] %||% "")
script_path <- if (nzchar(script_arg)) {
  normalizePath(script_arg, winslash = "/", mustWork = TRUE)
} else {
  normalizePath("papers/gut-dcst/scripts/build_ibd_reviewer_sensitivity_assets.R", winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_path)
paper_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
current_projects_root <- normalizePath(file.path(paper_root, "..", "..", ".."), winslash = "/", mustWork = TRUE)
gut_root <- file.path(current_projects_root, "gut_microbiome")
tables_dir <- file.path(paper_root, "assets", "tables")

run_dir <- file.path(gut_root, "outputs", "dcst_analysis", "runs", "2026-04-11-absorb-depthscan-adaptive")
validation_dir <- file.path(gut_root, "outputs", "dcst_validation")
agp_metadata_path <- file.path(gut_root, "data", "prime_gut_project_sample_metadata_2026-03-24.csv.gz")

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

fmt_num <- function(x, digits = 3) {
  if (is.na(x) || !is.finite(x)) {
    return("NA")
  }
  if (abs(x) < 0.001 && x != 0) {
    return(sprintf("%.2e", x))
  }
  sprintf(paste0("%.", digits, "g"), x)
}

clean_status <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y == "" | tolower(y) == "nan"] <- NA_character_
  y <- ifelse(tolower(y) == "year", "Year", y)
  y <- ifelse(tolower(y) == "week", "Week", y)
  y <- ifelse(tolower(y) == "month", "Month", y)
  y
}

safe_fisher <- function(a, b, c, d) {
  tab <- matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE)
  if (any(rowSums(tab) == 0L) || any(colSums(tab) == 0L)) {
    return(list(or = NA_real_, p = NA_real_))
  }
  fit <- fisher.test(tab)
  list(or = unname(fit$estimate), p = unname(fit$p.value))
}

safe_glm <- function(df) {
  fit <- tryCatch(
    suppressWarnings(glm(outcome ~ label_present + Host_Age + Host_Sex + Host_BMI, data = df, family = binomial())),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, status = fit$message))
  }
  tab <- coef(summary(fit))
  if (!"label_present" %in% rownames(tab)) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, status = "coefficient unavailable"))
  }
  beta <- unname(tab["label_present", "Estimate"])
  se <- unname(tab["label_present", "Std. Error"])
  p <- unname(tab["label_present", "Pr(>|z|)"])
  if (!is.finite(beta) || !is.finite(se) || se <= 0) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, status = "non-finite estimate"))
  }
  list(or = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se), p = p, status = "ok")
}

safe_glm_antibiotic <- function(df) {
  fit <- tryCatch(
    suppressWarnings(glm(outcome ~ label_present + Host_Age + Host_Sex + Host_BMI + Antibiotic_Status_Group, data = df, family = binomial())),
    error = function(e) e
  )
  if (inherits(fit, "error")) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, status = fit$message))
  }
  tab <- coef(summary(fit))
  if (!"label_present" %in% rownames(tab)) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, status = "coefficient unavailable"))
  }
  beta <- unname(tab["label_present", "Estimate"])
  se <- unname(tab["label_present", "Std. Error"])
  p <- unname(tab["label_present", "Pr(>|z|)"])
  if (!is.finite(beta) || !is.finite(se) || se <= 0) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_, p = NA_real_, status = "non-finite estimate"))
  }
  list(or = exp(beta), lo = exp(beta - 1.96 * se), hi = exp(beta + 1.96 * se), p = p, status = "ok")
}

prep_covariates <- function(df) {
  age_levels <- c("baby", "child", "teen", "20-29", "30-39", "40-49", "50-59", "60-69", "above 70")
  out <- df
  out$outcome <- as.integer(out$IBD)
  out$Host_BMI <- suppressWarnings(as.numeric(out$Host_BMI))
  out$Host_Sex <- factor(as.character(out$Host_Sex), levels = c("Female", "Male"))
  out$Host_Age <- factor(as.character(out$Host_Age), levels = age_levels)
  out$Antibiotic_Status_Group <- factor(clean_status(out$Antibiotic_Status), levels = c("No", "Year", "6 months", "Month", "Week"))
  out
}

label_specs <- data.frame(
  label_short = c(
    "Bacteroides",
    "Bacteroides / Lachnospiraceae",
    "Bacteroides / Lachnospiraceae / Alistipes",
    "Bacteroides / Lachnospiraceae / Alistipes / Faecalibacterium",
    "Morganella",
    "Proteus"
  ),
  depth = c(1L, 2L, 3L, 4L, 1L, 1L),
  label = c(
    "Bacteroides",
    "Bacteroides__Lachnospiraceae",
    "Bacteroides__Lachnospiraceae__Alistipes",
    "Bacteroides__Lachnospiraceae__Alistipes__Faecalibacterium",
    "Morganella",
    "Proteus"
  ),
  stringsAsFactors = FALSE
)

build_antibiotic_sensitivity <- function() {
  assign <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(run_dir, "agp_absorb_assignments.tsv.gz"))), data.table = FALSE)
  meta <- fread(
    cmd = sprintf("gzip -dc %s", shQuote(agp_metadata_path)),
    select = c("Run", "BioProject", "Antibiotic_Status"),
    data.table = FALSE
  )
  meta <- meta[meta$BioProject == "PRJEB11419", c("Run", "Antibiotic_Status"), drop = FALSE]
  assign$Antibiotic_Status <- meta$Antibiotic_Status[match(assign$Run, meta$Run)]
  assign <- prep_covariates(assign)

  rows <- list()
  for (i in seq_len(nrow(label_specs))) {
    spec <- label_specs[i, ]
    label_col <- sprintf("depth%d_absorb", spec$depth)
    base <- assign
    base$label_present <- as.integer(base[[label_col]] == spec$label)

    complete <- base[complete.cases(base[, c("outcome", "label_present", "Host_Age", "Host_Sex", "Host_BMI")]), , drop = FALSE]
    known_abx <- complete[!is.na(complete$Antibiotic_Status_Group), , drop = FALSE]
    no_recent <- known_abx[!as.character(known_abx$Antibiotic_Status_Group) %in% c("Week", "Month", "6 months"), , drop = FALSE]

    fit_primary <- safe_glm(complete)
    fit_abx <- safe_glm_antibiotic(known_abx)
    fit_no_recent <- safe_glm(no_recent)

    for (analysis_name in c("age_sex_bmi_complete_case", "plus_antibiotic_status", "excluding_week_month_6mo_antibiotics")) {
      df_model <- switch(
        analysis_name,
        age_sex_bmi_complete_case = complete,
        plus_antibiotic_status = known_abx,
        excluding_week_month_6mo_antibiotics = no_recent
      )
      fit <- switch(
        analysis_name,
        age_sex_bmi_complete_case = fit_primary,
        plus_antibiotic_status = fit_abx,
        excluding_week_month_6mo_antibiotics = fit_no_recent
      )
      rows[[length(rows) + 1L]] <- data.frame(
        label_short = spec$label_short,
        depth = spec$depth,
        analysis = analysis_name,
        n = nrow(df_model),
        cases = sum(df_model$outcome == 1L),
        controls = sum(df_model$outcome == 0L),
        label_cases = sum(df_model$outcome == 1L & df_model$label_present == 1L),
        label_controls = sum(df_model$outcome == 0L & df_model$label_present == 1L),
        adjusted_or = fit$or,
        ci_low = fit$lo,
        ci_high = fit$hi,
        p_value = fit$p,
        fit_status = fit$status,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, q_value_selected_labels := p.adjust(p_value, method = "BH"), by = analysis]
  fwrite(out, file.path(tables_dir, "TABLE_S9_ibd_antibiotic_status_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")
  out
}

association_table <- function(labels, condition) {
  labs <- sort(unique(labels[!is.na(labels) & nzchar(labels)]))
  rows <- lapply(labs, function(lab) {
    present <- labels == lab
    present[is.na(present)] <- FALSE
    a <- sum(present & condition == 1L)
    b <- sum(present & condition == 0L)
    c <- sum(!present & condition == 1L)
    d <- sum(!present & condition == 0L)
    fit <- safe_fisher(a, b, c, d)
    data.frame(label = lab, n_label = a + b, n_cases = a, n_controls = b, or = fit$or, p_value = fit$p, stringsAsFactors = FALSE)
  })
  out <- rbindlist(rows, fill = TRUE)
  if (nrow(out)) {
    out$q_value <- p.adjust(out$p_value, method = "BH")
    out <- out[order(out$q_value, out$p_value, out$label), ]
  }
  out
}

derive_condition <- function(df, comparison) {
  text_cols <- intersect(c("final_case_control_label", "phenotype_group", "phenotype_hint", "case_control_hint"), names(df))
  combined <- apply(df[, text_cols, drop = FALSE], 1, function(x) paste(tolower(as.character(x)), collapse = " "))
  negative <- grepl("healthy", combined)
  if (comparison == "Crohn_vs_UC") {
    positive <- grepl("crohn", combined)
  } else {
    positive <- !negative
  }
  as.integer(positive)
}

read_best_signal <- function(root, cohort, comparison) {
  best <- NULL
  for (depth in 1:4) {
    path <- file.path(root, cohort, sprintf("%s__%s_depth%d_results.csv", cohort, comparison, depth))
    if (!file.exists(path)) next
    df <- read.csv(path, check.names = FALSE)
    if (!nrow(df)) next
    df <- df[order(df$q_value, df$p_value, df$DCST), , drop = FALSE]
    row <- df[1, , drop = FALSE]
    row$depth <- depth
    if (is.null(best) || row$q_value[[1]] < best$q_value[[1]]) {
      best <- row
    }
  }
  best
}

build_external_subject_sensitivity <- function() {
  specs <- data.frame(
    cohort = c("halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "hmp2", "gevers_2014", "prjeb84421"),
    cohort_short = c("Halfvarson IBD", "Halfvarson Crohn", "Halfvarson UC", "HMP2 IBD", "Gevers Crohn", "PRJEB84421 OFG"),
    comparison = c("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy", "IBD_vs_Healthy", "Crohn_vs_Healthy", "OFG_vs_Healthy"),
    stringsAsFactors = FALSE
  )
  modes <- data.frame(
    mode = c("rebuilt_cohort", "agp_label_transfer"),
    root = c(file.path(validation_dir, "absorb"), file.path(validation_dir, "frozen_agp")),
    label_template = c("depth%d_absorb", "depth%d_frozen_agp_absorb"),
    stringsAsFactors = FALSE
  )
  rows <- list()
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    for (j in seq_len(nrow(modes))) {
      mode <- modes[j, ]
      best <- read_best_signal(mode$root, spec$cohort, spec$comparison)
      if (is.null(best)) next
      assign_path <- file.path(mode$root, spec$cohort, sprintf("%s__%s_dcst_assignments.csv", spec$cohort, spec$comparison))
      if (!file.exists(assign_path)) next
      assign <- read.csv(assign_path, check.names = FALSE)
      if (!"subject_id" %in% names(assign)) {
        assign$subject_id <- assign$run_accession
      }
      assign$subject_id <- as.character(assign$subject_id)
      assign$condition <- derive_condition(assign, spec$comparison)
      label_col <- sprintf(mode$label_template, best$depth[[1]])
      assign <- assign[order(assign$subject_id, assign$run_accession), , drop = FALSE]
      first <- assign[!duplicated(assign$subject_id), , drop = FALSE]
      labels <- as.character(first[[label_col]])
      mapped <- !is.na(labels) & nzchar(labels)
      assoc <- association_table(labels[mapped], first$condition[mapped])
      match_row <- assoc[assoc$label == as.character(best$DCST[[1]]), , drop = FALSE]
      if (!nrow(match_row)) {
        match_row <- data.frame(label = as.character(best$DCST[[1]]), n_label = NA_integer_, n_cases = NA_integer_, n_controls = NA_integer_, or = NA_real_, p_value = NA_real_, q_value = NA_real_)
      }
      rows[[length(rows) + 1L]] <- data.frame(
        cohort = spec$cohort_short,
        comparison = spec$comparison,
        mode = mode$mode,
        depth = best$depth[[1]],
        sample_best_label = as.character(best$DCST[[1]]),
        sample_or = best$OR[[1]],
        sample_q = best$q_value[[1]],
        n_samples = nrow(assign),
        n_subjects = length(unique(assign$subject_id)),
        repeated_subjects = sum(table(assign$subject_id) > 1L),
        first_subject_mapped_n = sum(mapped),
        first_subject_label_n = match_row$n_label[[1]],
        first_subject_label_cases = match_row$n_cases[[1]],
        first_subject_label_controls = match_row$n_controls[[1]],
        first_subject_or = match_row$or[[1]],
        first_subject_q = match_row$q_value[[1]],
        stringsAsFactors = FALSE
      )
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S10_external_subject_level_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")
  out
}

build_mapping_balance <- function() {
  specs <- data.frame(
    cohort = c("halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "hmp2", "gevers_2014", "prjeb84421", "jacobs_2023_ibs_250bp"),
    cohort_short = c("Halfvarson IBD", "Halfvarson Crohn", "Halfvarson UC", "HMP2 IBD", "Gevers Crohn", "PRJEB84421 OFG", "Jacobs IBS"),
    comparison = c("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy", "IBD_vs_Healthy", "Crohn_vs_Healthy", "OFG_vs_Healthy", "IBS_vs_Healthy"),
    stringsAsFactors = FALSE
  )
  rows <- list()
  for (i in seq_len(nrow(specs))) {
    spec <- specs[i, ]
    path <- file.path(validation_dir, "frozen_agp", spec$cohort, sprintf("%s__%s_mapping_summary.csv", spec$cohort, spec$comparison))
    if (!file.exists(path)) next
    df <- read.csv(path, check.names = FALSE)
    df$cohort <- spec$cohort_short
    df$comparison <- spec$comparison
    rows[[length(rows) + 1L]] <- df
  }
  out <- rbindlist(rows, fill = TRUE)
  keep <- c("cohort", "comparison", "depth", "mapped_total", "mapped_cases", "mapped_controls", "mapping_rate_total", "mapping_rate_cases", "mapping_rate_controls", "mapped_vs_unmapped_p")
  out <- out[, keep, with = FALSE]
  fwrite(out, file.path(tables_dir, "TABLE_S11_external_label_transfer_mapping_balance.tsv"), sep = "\t", quote = FALSE, na = "NA")
  out
}

build_absorb_policy_sensitivity <- function() {
  obj <- readRDS(file.path(run_dir, "agp_absorb_object.rds"))
  assign <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(run_dir, "agp_absorb_assignments.tsv.gz"))), data.table = FALSE)
  ibd_by_run <- setNames(as.integer(assign$IBD), assign$Run)
  rows <- list()
  for (i in seq_len(nrow(label_specs))) {
    spec <- label_specs[i, ]
    for (view in c("absorb", "rare_policy")) {
      labels <- if (view == "absorb") obj$cst.levels.absorb[[spec$depth]] else obj$cst.levels.rare[[spec$depth]]
      condition <- ibd_by_run[names(labels)]
      present <- as.character(labels) == spec$label
      present[is.na(present)] <- FALSE
      a <- sum(present & condition == 1L, na.rm = TRUE)
      b <- sum(present & condition == 0L, na.rm = TRUE)
      c <- sum(!present & condition == 1L, na.rm = TRUE)
      d <- sum(!present & condition == 0L, na.rm = TRUE)
      fit <- safe_fisher(a, b, c, d)
      rows[[length(rows) + 1L]] <- data.frame(
        label_short = spec$label_short,
        depth = spec$depth,
        view = view,
        n_label = a + b,
        label_cases = a,
        label_controls = b,
        fisher_or = fit$or,
        p_value = fit$p,
        stringsAsFactors = FALSE
      )
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  out[, q_value_selected_labels := p.adjust(p_value, method = "BH"), by = view]
  fwrite(out, file.path(tables_dir, "TABLE_S12_absorb_policy_selected_ibd_robustness.tsv"), sep = "\t", quote = FALSE, na = "NA")
  out
}

parse_rank <- function(taxonomy) {
  ranks <- strsplit(as.character(taxonomy), ";", fixed = TRUE)[[1L]]
  pick <- function(prefix) {
    hit <- ranks[startsWith(ranks, prefix)]
    if (!length(hit)) return("")
    out <- sub(paste0("^", prefix), "", hit[[1L]])
    out <- trimws(out)
    if (!nzchar(out) || grepl("^_+$", out)) "" else out
  }
  if (nzchar(pick("s__"))) return("species")
  if (nzchar(pick("g__"))) return("genus")
  if (nzchar(pick("f__"))) return("family")
  if (nzchar(pick("o__"))) return("order")
  if (nzchar(pick("c__"))) return("class")
  if (nzchar(pick("p__"))) return("phylum")
  if (nzchar(pick("d__"))) return("domain")
  "unresolved"
}

build_taxonomic_rank_table <- function() {
  ref <- fread(file.path(run_dir, "taxon_reference.tsv"), data.table = FALSE)
  components <- unique(unlist(strsplit(label_specs$label, "__", fixed = TRUE)))
  rows <- lapply(components, function(component) {
    hit <- ref[ref$feature_label == component, , drop = FALSE]
    if (!nrow(hit)) {
      return(data.frame(component = component, taxonomic_rank = "not_found", taxonomy = "", stringsAsFactors = FALSE))
    }
    data.frame(component = component, taxonomic_rank = parse_rank(hit$taxonomy[[1]]), taxonomy = hit$taxonomy[[1]], stringsAsFactors = FALSE)
  })
  out <- rbindlist(rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S13_selected_ibd_label_taxonomic_ranks.tsv"), sep = "\t", quote = FALSE, na = "NA")
  out
}

write_summary <- function(antibiotic, subject, mapping, policy, ranks) {
  lines <- c(
    "# IBD Reviewer Sensitivity Asset Summary",
    "",
    "Generated by `build_ibd_reviewer_sensitivity_assets.R`.",
    "",
    "Generated manuscript-facing tables:",
    "- `TABLE_S9_ibd_antibiotic_status_sensitivity.tsv`",
    "- `TABLE_S10_external_subject_level_sensitivity.tsv`",
    "- `TABLE_S11_external_label_transfer_mapping_balance.tsv`",
    "- `TABLE_S12_absorb_policy_selected_ibd_robustness.tsv`",
    "- `TABLE_S13_selected_ibd_label_taxonomic_ranks.tsv`",
    "",
    sprintf("Antibiotic-status sensitivity rows: %d", nrow(antibiotic)),
    sprintf("External subject-level sensitivity rows: %d", nrow(subject)),
    sprintf("External mapping-balance rows: %d", nrow(mapping)),
    sprintf("Absorb-policy sensitivity rows: %d", nrow(policy)),
    sprintf("Selected component taxonomic-rank rows: %d", nrow(ranks))
  )
  writeLines(lines, file.path(tables_dir, "TABLE_S9_S13_reviewer_sensitivity_summary.md"))
}

antibiotic <- build_antibiotic_sensitivity()
subject <- build_external_subject_sensitivity()
mapping <- build_mapping_balance()
policy <- build_absorb_policy_sensitivity()
ranks <- build_taxonomic_rank_table()
write_summary(antibiotic, subject, mapping, policy, ranks)

message("Wrote reviewer sensitivity assets to ", tables_dir)
