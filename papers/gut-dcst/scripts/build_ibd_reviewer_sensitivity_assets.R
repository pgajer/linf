#!/usr/bin/env Rscript

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
  normalizePath("papers/gut-dcst/scripts/build_ibd_reviewer_sensitivity_assets.R", winslash = "/", mustWork = TRUE)
}
script_dir <- dirname(script_path)
paper_root <- normalizePath(file.path(script_dir, ".."), winslash = "/", mustWork = TRUE)
current_projects_root <- normalizePath(file.path(paper_root, "..", "..", ".."), winslash = "/", mustWork = TRUE)
gut_root <- file.path(current_projects_root, "gut_microbiome")
tables_dir <- file.path(paper_root, "assets", "tables")

run_dir <- file.path(gut_root, "outputs", "dcst_analysis", "runs", "2026-04-26-agp-silva-local-qza-absorb-depth4-n0_50_25_25_25")
validation_dir <- file.path(gut_root, "outputs", "dcst_validation")
transfer_root <- file.path(validation_dir, "frozen_agp_silva_local_qza_2026-04-26_n0_50_25_25_25")
agp_metadata_path <- file.path(
  gut_root, "outputs", "agp_silva_taxonomy", "2026-04-26-redbiom-md5-sequence-map",
  "agp_all_good_depthscan_metadata.tsv.gz"
)
agp_counts_path <- file.path(
  gut_root, "outputs", "agp_silva_taxonomy", "2026-04-26-redbiom-md5-sequence-map",
  "agp_silva_collapsed_counts.tsv.gz"
)

dir.create(tables_dir, recursive = TRUE, showWarnings = FALSE)

bootstrap_iterations <- max(1L, suppressWarnings(as.integer(Sys.getenv("DCST_BOOTSTRAP_B", "1000"))))
bootstrap_seed <- suppressWarnings(as.integer(Sys.getenv("DCST_BOOTSTRAP_SEED", "20260425")))
if (!is.finite(bootstrap_seed)) {
  bootstrap_seed <- 20260425L
}

tex_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("\\\\", "\\\\textbackslash{}", x)
  x <- gsub("([%&_#{}$])", "\\\\\\1", x, perl = TRUE)
  x
}

fmt_prob <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("NA")
  }
  if (x < 0.001) {
    return(formatC(x, format = "e", digits = 1))
  }
  formatC(x, format = "f", digits = 3)
}

fmt_or <- function(x) {
  if (is.na(x)) {
    return("NA")
  }
  if (is.infinite(x)) {
    return(ifelse(x > 0, "inf", "-inf"))
  }
  trimws(formatC(x, format = "fg", digits = 3))
}

fmt_pct <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("NA")
  }
  formatC(x, format = "f", digits = 1)
}

fmt_int <- function(x) {
  if (is.na(x) || !is.finite(x)) {
    return("NA")
  }
  formatC(x, format = "d", big.mark = ",")
}

stable_seed_from_text <- function(text) {
  ints <- utf8ToInt(as.character(text))
  if (!length(ints)) {
    return(1L)
  }
  weights <- seq_along(ints)
  as.integer((sum(ints * weights) %% 1000000000L) + 1L)
}

bayes_or_ci <- function(a, b, c, d, seed_key, draws = 20000L) {
  set.seed(stable_seed_from_text(seed_key))
  p_case <- rbeta(draws, a + 0.5, c + 0.5)
  p_control <- rbeta(draws, b + 0.5, d + 0.5)
  odds_case <- p_case / (1 - p_case)
  odds_control <- p_control / (1 - p_control)
  or_draws <- odds_case / odds_control
  finite <- is.finite(or_draws)
  if (!any(finite)) {
    return(list(or = NA_real_, lo = NA_real_, hi = NA_real_))
  }
  keep <- or_draws[finite]
  list(
    or = median(keep),
    lo = as.numeric(quantile(keep, probs = 0.025)),
    hi = as.numeric(quantile(keep, probs = 0.975))
  )
}

clean_status <- function(x) {
  y <- trimws(as.character(x))
  y[is.na(y) | y == "" | tolower(y) == "nan"] <- NA_character_
  y <- ifelse(tolower(y) == "year", "Year", y)
  y <- ifelse(tolower(y) == "week", "Week", y)
  y <- ifelse(tolower(y) == "month", "Month", y)
  y
}

normalize_age_group <- function(x) {
  y <- trimws(tolower(as.character(x)))
  y[y %in% c("", "na", "nan", "not provided", "not applicable", "not collected", "labcontrol test")] <- NA_character_
  y <- ifelse(y == "20s", "20-29", y)
  y <- ifelse(y == "30s", "30-39", y)
  y <- ifelse(y == "40s", "40-49", y)
  y <- ifelse(y == "50s", "50-59", y)
  y <- ifelse(y == "60s", "60-69", y)
  y <- ifelse(y == "70+", "above 70", y)
  y
}

normalize_sex_group <- function(x) {
  y <- trimws(tolower(as.character(x)))
  y[y %in% c("", "na", "nan", "not provided", "not applicable", "not collected", "other")] <- NA_character_
  y <- ifelse(y == "female", "Female", y)
  y <- ifelse(y == "male", "Male", y)
  y
}

normalize_binary_flag <- function(x) {
  y <- trimws(tolower(as.character(x)))
  y[is.na(y) | y == "" | y %in% c("nan", "na", "missing", "unknown", "not applicable", "not collected")] <- NA_character_
  y[!(y %in% c("true", "false") | is.na(y))] <- NA_character_
  y
}

derive_agp_antibiotic_status <- function(x) {
  y <- trimws(tolower(as.character(x)))
  y[is.na(y) | y == "" | y %in% c("nan", "na", "missing", "unknown")] <- NA_character_
  out <- rep(NA_character_, length(y))
  out[grepl("have not taken antibiotics in the past year|no antibiotics in the past year", y)] <- "No"
  out[grepl("past week|\\bweek\\b", y)] <- "Week"
  out[grepl("past month|\\bmonth\\b", y)] <- "Month"
  out[grepl("6 months|six months", y)] <- "6 months"
  out[grepl("\\byear\\b", y) & is.na(out)] <- "Year"
  out[grepl("not taken antibiotics", y) & is.na(out)] <- "No"
  out
}

safe_fisher <- function(a, b, c, d) {
  tab <- matrix(c(a, b, c, d), nrow = 2L, byrow = TRUE)
  if (any(rowSums(tab) == 0L) || any(colSums(tab) == 0L)) {
    return(list(or = NA_real_, p = NA_real_, ci_low = NA_real_, ci_high = NA_real_))
  }
  fit <- fisher.test(tab)
  list(
    or = unname(fit$estimate),
    p = unname(fit$p.value),
    ci_low = unname(fit$conf.int[[1L]]),
    ci_high = unname(fit$conf.int[[2L]])
  )
}

haldane_or <- function(a, b, c, d) {
  ((a + 0.5) * (d + 0.5)) / ((b + 0.5) * (c + 0.5))
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
  out$Host_Sex <- factor(normalize_sex_group(out$Host_Sex), levels = c("Female", "Male"))
  out$Host_Age <- factor(normalize_age_group(out$Host_Age), levels = age_levels)
  if (!"Antibiotic_Status" %in% names(out)) {
    out$Antibiotic_Status <- NA_character_
  }
  out$Antibiotic_Status_Group <- factor(clean_status(out$Antibiotic_Status), levels = c("No", "Year", "6 months", "Month", "Week"))
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

summarize_taxonomy <- function(taxonomy) {
  parts <- strsplit(as.character(taxonomy), ";", fixed = TRUE)[[1L]]
  keep <- trimws(sub("^[a-z]__", "", parts))
  keep <- keep[nzchar(keep) & keep != "_"]
  if (!length(keep)) {
    return("")
  }
  paste(keep[seq_len(min(length(keep), 5L))], collapse = "; ")
}

taxon_reference <- fread(file.path(run_dir, "taxon_reference.tsv"), data.table = FALSE)
taxonomy_lookup <- setNames(taxon_reference$taxonomy, taxon_reference$feature_label)
rank_lookup <- setNames(vapply(taxon_reference$taxonomy, parse_rank, character(1L)), taxon_reference$feature_label)

format_component_tex <- function(component) {
  component <- as.character(component)
  if (!nzchar(component) || is.na(component)) {
    return("NA")
  }
  if (component == "RARE_DOMINANT") {
    return("\\emph{RARE\\_DOMINANT}")
  }
  rank <- if (component %in% names(rank_lookup)) rank_lookup[[component]] else NA_character_
  if (is.na(rank) || !nzchar(rank)) {
    rank <- if (grepl("aceae$", component)) {
      "family"
    } else if (grepl(" ", component, fixed = TRUE)) {
      "species"
    } else {
      "genus"
    }
  }
  escaped <- tex_escape(component)
  if (rank == "family") {
    return(escaped)
  }
  if (rank == "species" || grepl(" ", component, fixed = TRUE)) {
    return(sprintf("\\emph{%s}", escaped))
  }
  sprintf("\\taxon{%s}", escaped)
}

format_label_tex <- function(label) {
  label <- as.character(label)
  if (is.na(label) || !nzchar(label)) {
    return("NA")
  }
  parts <- strsplit(label, "__", fixed = TRUE)[[1L]]
  paste(vapply(parts, format_component_tex, character(1L)), collapse = "\\dcstsep ")
}

comparison_display <- function(comparison) {
  switch(
    comparison,
    IBD_vs_Healthy = "IBD vs healthy",
    Crohn_vs_Healthy = "Crohn vs healthy",
    UC_vs_Healthy = "UC vs healthy",
    OFG_vs_Healthy = "OFG vs healthy",
    IBS_vs_Healthy = "IBS vs healthy",
    Crohn_vs_UC = "Crohn vs UC",
    tex_escape(gsub("_", " ", comparison, fixed = TRUE))
  )
}

interpret_overview <- function(mode, comparison, q_value, coverage = NA_real_) {
  if (is.na(q_value)) {
    return("No evaluable signal.")
  }
  if (mode == "rebuilt") {
    if (q_value < 0.05) {
      if (comparison == "UC_vs_Healthy") return("Corrected rebuilt UC signal detected.")
      if (comparison == "Crohn_vs_Healthy") return("Corrected rebuilt Crohn signal detected.")
      if (comparison == "IBD_vs_Healthy") return("Corrected rebuilt pooled-IBD signal detected.")
      return("Corrected rebuilt signal detected.")
    }
    if (q_value < 0.25) {
      return("Directional only in rebuilt mode.")
    }
    return("No corrected rebuilt signal.")
  }
  if (q_value < 0.05) {
    if (comparison == "UC_vs_Healthy") return("Corrected AGP-derived UC transfer detected.")
    if (comparison == "Crohn_vs_Healthy") return("Corrected AGP-derived Crohn transfer detected.")
    if (comparison == "IBD_vs_Healthy") return("Corrected AGP-derived pooled-IBD transfer detected.")
    return("Corrected AGP-derived transfer detected.")
  }
  if (is.finite(coverage) && coverage >= 75) {
    return("Samples map into the hierarchy, but no corrected transferred signal.")
  }
  "Limited mapping constrains transferred signal detection."
}

interpret_mapping_balance <- function(case_rate, control_rate, p_value, depth) {
  if (!is.finite(case_rate) || !is.finite(control_rate) || !is.finite(p_value)) {
    return("Unavailable")
  }
  diff <- case_rate - control_rate
  if (p_value < 0.05 && diff > 0) {
    if (depth >= 4) return("Depth-4 case skew")
    return("Cases higher")
  }
  if (p_value < 0.05 && diff < 0) {
    return("Controls higher")
  }
  if (abs(diff) <= 5) {
    return("Approx. balanced")
  }
  "No strong imbalance"
}

followup_score <- function(dt) {
  score <- rep(-Inf, nrow(dt))
  q_ok <- is.finite(dt$adj_freq_q_global) & !is.na(dt$adj_freq_q_global)
  score[q_ok] <- 100 + pmax(0, -log10(pmax(dt$adj_freq_q_global[q_ok], 1e-12)))
  bayes_ok <- is.finite(dt$bayes_adj_or) & !is.na(dt$bayes_adj_or) & dt$bayes_adj_or > 0
  score[!q_ok & bayes_ok] <- abs(log(dt$bayes_adj_or[!q_ok & bayes_ok]))
  score
}

build_dynamic_label_specs <- function() {
  followup_path <- file.path(run_dir, "label_followup", "ibd", "ibd_label_level_results.tsv")
  dt <- fread(followup_path)
  dt[, score := followup_score(dt)]
  dt <- dt[is.finite(score)]
  setorder(dt, depth, -score, -n_label_all)

  chosen <- vector("list", 0L)
  chosen_labels <- character()

  for (d in 1:4) {
    sub <- dt[depth == d]
    if (!nrow(sub)) {
      next
    }
    pick_idx <- which(!sub$label %in% chosen_labels)[1]
    if (!length(pick_idx) || is.na(pick_idx)) {
      pick_idx <- 1L
    }
    row <- sub[pick_idx]
    chosen[[length(chosen) + 1L]] <- row
    chosen_labels <- c(chosen_labels, row$label)
  }

  extras <- copy(dt)
  setorder(extras, -score, -n_label_all)
  for (i in seq_len(nrow(extras))) {
    row <- extras[i]
    if (row$label %in% chosen_labels) {
      next
    }
    chosen[[length(chosen) + 1L]] <- row
    chosen_labels <- c(chosen_labels, row$label)
    if (length(chosen) >= 6L) {
      break
    }
  }

  out <- rbindlist(chosen, fill = TRUE)
  setorder(out, depth, -score, -n_label_all)
  data.frame(
    label_short = gsub("__", " / ", out$label, fixed = TRUE),
    depth = as.integer(out$depth),
    label = out$label,
    stringsAsFactors = FALSE
  )
}

label_specs <- build_dynamic_label_specs()

build_cleaner_control_sensitivity <- function() {
  assign <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(run_dir, "agp_absorb_assignments.tsv.gz"))), data.table = FALSE)
  assign <- prep_covariates(assign)
  cleaner_mask <- with(assign, IBD == 1L | (IBD == 0L & IBS == 0L & Autoimmune == 0L & Acid_reflux == 0L & Seasonal_allergies == 0L))
  cleaner <- assign[cleaner_mask, , drop = FALSE]

  rows <- list()
  for (i in seq_len(nrow(label_specs))) {
    spec <- label_specs[i, ]
    label_col <- sprintf("depth%d_absorb", spec$depth)

    primary <- assign
    primary$label_present <- as.integer(primary[[label_col]] == spec$label)
    primary_complete <- primary[complete.cases(primary[, c("outcome", "label_present", "Host_Age", "Host_Sex", "Host_BMI")]), , drop = FALSE]
    primary_fit <- safe_glm(primary_complete)

    clean <- cleaner
    clean$label_present <- as.integer(clean[[label_col]] == spec$label)
    clean_complete <- clean[complete.cases(clean[, c("outcome", "label_present", "Host_Age", "Host_Sex", "Host_BMI")]), , drop = FALSE]
    clean_fit <- safe_glm(clean_complete)

    rows[[length(rows) + 1L]] <- data.frame(
      label_short = spec$label_short,
      depth = spec$depth,
      label = spec$label,
      primary_cases = sum(primary$outcome == 1L & primary$label_present == 1L, na.rm = TRUE),
      primary_controls = sum(primary$outcome == 0L & primary$label_present == 1L, na.rm = TRUE),
      primary_adj_or = primary_fit$or,
      primary_p = primary_fit$p,
      primary_status = primary_fit$status,
      cleaner_cases = sum(clean$outcome == 1L & clean$label_present == 1L, na.rm = TRUE),
      cleaner_controls = sum(clean$outcome == 0L & clean$label_present == 1L, na.rm = TRUE),
      cleaner_adj_or = clean_fit$or,
      cleaner_p = clean_fit$p,
      cleaner_status = clean_fit$status,
      stringsAsFactors = FALSE
    )
  }

  out <- rbindlist(rows, fill = TRUE)
  out[, primary_q := p.adjust(primary_p, method = "BH")]
  out[, cleaner_q := p.adjust(cleaner_p, method = "BH")]
  out[, interpretation := ifelse(
    is.finite(primary_q) & primary_q < 0.05 & is.finite(cleaner_q) & cleaner_q < 0.05,
    "Signal retained after cleaner-control restriction.",
    ifelse(
      is.finite(primary_adj_or) & is.finite(cleaner_adj_or) &
        ((primary_adj_or > 1 & cleaner_adj_or > 1) | (primary_adj_or < 1 & cleaner_adj_or < 1)),
      "Direction retained, but support attenuates after cleaner-control restriction.",
      "Cleaner-control restriction does not support the same label pattern."
    )
  )]

  fwrite(out, file.path(tables_dir, "TABLE_S6_agp_ibd_cleaner_control_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")

  render_rows <- vapply(seq_len(nrow(out)), function(i) {
    row <- out[i, ]
    sprintf(
      "%s & %d & %s/%s; adjusted OR %s, q = %s. & %s/%s; adjusted OR %s, q = %s. & %s \\\\",
      format_label_tex(row$label),
      row$depth,
      fmt_int(row$primary_cases),
      fmt_int(row$primary_controls),
      fmt_or(row$primary_adj_or),
      fmt_prob(row$primary_q),
      fmt_int(row$cleaner_cases),
      fmt_int(row$cleaner_controls),
      fmt_or(row$cleaner_adj_or),
      fmt_prob(row$cleaner_q),
      tex_escape(row$interpretation)
    )
  }, character(1L))
  render_table(
    file.path(tables_dir, "TABLE_S6_agp_ibd_cleaner_control_sensitivity.tex"),
    "Y r Y Y Y",
    c("Label", "Depth", "Primary analysis", "Cleaner-control sensitivity", "Interpretation"),
    render_rows,
    size = "\\footnotesize"
  )
  out
}

build_retained_feature_profiles <- function() {
  assign <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(run_dir, "agp_absorb_assignments.tsv.gz"))), data.table = FALSE)
  counts <- fread(cmd = sprintf("gzip -dc %s", shQuote(agp_counts_path)), data.table = FALSE)
  counts <- as.data.frame(counts)
  rownames(counts) <- counts$Run
  counts$Run <- NULL

  rows <- list()
  for (i in seq_len(nrow(label_specs))) {
    spec <- label_specs[i, ]
    label_col <- sprintf("depth%d_absorb", spec$depth)
    mask <- assign[[label_col]] == spec$label
    sub_assign <- assign[mask, , drop = FALSE]
    if (!nrow(sub_assign)) {
      next
    }
    components <- strsplit(spec$label, "__", fixed = TRUE)[[1L]]
    keep_components <- components[components %in% colnames(counts)]
    if (!length(keep_components)) {
      next
    }
    rel <- counts[sub_assign$Run, keep_components, drop = FALSE]
    rel <- rel / rowSums(rel)
    row_all <- vapply(keep_components, function(comp) median(rel[[comp]], na.rm = TRUE) * 100, numeric(1L))
    ibd_runs <- sub_assign$Run[sub_assign$IBD == 1L]
    rel_ibd <- if (length(ibd_runs)) rel[ibd_runs, , drop = FALSE] else rel[0, , drop = FALSE]
    row_ibd <- vapply(keep_components, function(comp) {
      if (!nrow(rel_ibd)) return(NA_real_)
      median(rel_ibd[[comp]], na.rm = TRUE) * 100
    }, numeric(1L))

    rows[[length(rows) + 1L]] <- data.frame(
      label_short = spec$label_short,
      depth = spec$depth,
      label = spec$label,
      n_ibd = sum(sub_assign$IBD == 1L, na.rm = TRUE),
      n_control = sum(sub_assign$IBD == 0L, na.rm = TRUE),
      ordered_components = paste(vapply(keep_components, format_component_tex, character(1L)), collapse = "; "),
      median_all_pct = paste(vapply(row_all, fmt_pct, character(1L)), collapse = "; "),
      median_ibd_pct = paste(vapply(row_ibd, fmt_pct, character(1L)), collapse = "; "),
      stringsAsFactors = FALSE
    )
  }

  out <- rbindlist(rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S7_ibd_dominance_lineage_retained_feature_profiles.tsv"), sep = "\t", quote = FALSE, na = "NA")

  render_rows <- vapply(seq_len(nrow(out)), function(i) {
    row <- out[i, ]
    sprintf(
      "%s & %d & %s/%s & %s & %s & %s \\\\",
      format_label_tex(row$label),
      row$depth,
      fmt_int(row$n_ibd),
      fmt_int(row$n_control),
      row$ordered_components,
      tex_escape(row$median_all_pct),
      tex_escape(row$median_ibd_pct)
    )
  }, character(1L))
  render_table(
    file.path(tables_dir, "TABLE_S7_ibd_dominance_lineage_retained_feature_profiles.tex"),
    "Y r r Y Y Y",
    c("Label", "Depth", "IBD/control", "Ordered retained features", "Median all (\\%)", "Median IBD (\\%)"),
    render_rows,
    size = "\\scriptsize"
  )
  out
}

overview_specs <- data.frame(
  cohort = c("halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "hmp2", "gevers_2014", "prjeb84421", "jacobs_2023_ibs_250bp"),
  cohort_short = c("Halfvarson 2017", "Halfvarson 2017", "Halfvarson 2017", "HMP2 / IBDMDB", "Gevers 2014", "PRJEB84421", "Jacobs 2023 IBS 250bp"),
  comparison = c("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy", "IBD_vs_Healthy", "Crohn_vs_Healthy", "OFG_vs_Healthy", "IBS_vs_Healthy"),
  context_prefix = c(
    "Stool, longitudinal IBD cohort",
    "Stool, longitudinal IBD cohort",
    "Stool, longitudinal IBD cohort",
    "Mucosal gut-site IBDMDB samples",
    "Pediatric stool Crohn disease comparison",
    "Stool-based pediatric-onset orofacial granulomatosis comparison",
    "250bp IBS validation subset"
  ),
  stringsAsFactors = FALSE
)

focus_specs <- data.frame(
  cohort = c("halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "hmp2", "hmp2", "hmp2", "gevers_2014"),
  cohort_short = c("Halfvarson 2017", "Halfvarson 2017", "Halfvarson 2017", "HMP2 / IBDMDB", "HMP2 / IBDMDB", "HMP2 / IBDMDB", "Gevers 2014"),
  comparison = c("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy", "IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy", "Crohn_vs_Healthy"),
  stringsAsFactors = FALSE
)

subject_diag_specs <- data.frame(
  cohort = c("halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "hmp2", "gevers_2014", "prjeb84421"),
  cohort_short = c("Halfvarson 2017", "Halfvarson 2017", "Halfvarson 2017", "HMP2 / IBDMDB", "Gevers 2014", "PRJEB84421"),
  comparison = c("IBD_vs_Healthy", "Crohn_vs_Healthy", "UC_vs_Healthy", "IBD_vs_Healthy", "Crohn_vs_Healthy", "OFG_vs_Healthy"),
  stringsAsFactors = FALSE
)

mapping_diag_specs <- data.frame(
  cohort = c("halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "halfvarson_2017", "hmp2", "hmp2", "gevers_2014", "prjeb84421", "jacobs_2023_ibs_250bp"),
  cohort_short = c("Halfvarson 2017", "Halfvarson 2017", "Halfvarson 2017", "Halfvarson 2017", "HMP2 / IBDMDB", "HMP2 / IBDMDB", "Gevers 2014", "PRJEB84421", "Jacobs 2023 IBS 250bp"),
  comparison = c("IBD_vs_Healthy", "IBD_vs_Healthy", "UC_vs_Healthy", "UC_vs_Healthy", "IBD_vs_Healthy", "IBD_vs_Healthy", "Crohn_vs_Healthy", "OFG_vs_Healthy", "IBS_vs_Healthy"),
  depth = c(2L, 4L, 2L, 4L, 2L, 4L, 2L, 2L, 2L),
  stringsAsFactors = FALSE
)

load_assignment_csv <- function(mode_root, cohort, comparison) {
  path <- file.path(mode_root, cohort, sprintf("%s__%s_dcst_assignments.csv", cohort, comparison))
  if (!file.exists(path)) {
    return(NULL)
  }
  out <- fread(path, data.table = FALSE)
  if (!"subject_id" %in% names(out)) {
    out$subject_id <- out$run_accession
  }
  out$subject_id <- as.character(out$subject_id)
  out
}

association_table <- function(labels, condition) {
  mapped <- !is.na(labels) & nzchar(labels)
  labels <- as.character(labels[mapped])
  condition <- as.integer(condition[mapped])
  labs <- sort(unique(labels))
  if (!length(labs)) {
    return(data.frame(label = character(), n_label = integer(), n_cases = integer(), n_controls = integer(), or = numeric(), p_value = numeric(), q_value = numeric(), stringsAsFactors = FALSE))
  }
  total_cases <- sum(condition == 1L)
  total_controls <- sum(condition == 0L)
  case_idx <- match(labels[condition == 1L], labs)
  ctrl_idx <- match(labels[condition == 0L], labs)
  counts_cases <- tabulate(case_idx, nbins = length(labs))
  counts_controls <- tabulate(ctrl_idx, nbins = length(labs))
  pvals <- numeric(length(labs))
  ors <- numeric(length(labs))
  ci_lows <- numeric(length(labs))
  ci_highs <- numeric(length(labs))
  for (k in seq_along(labs)) {
    a <- counts_cases[[k]]
    b <- counts_controls[[k]]
    c <- total_cases - a
    d <- total_controls - b
    fit <- safe_fisher(a, b, c, d)
    pvals[[k]] <- fit$p
    ors[[k]] <- fit$or
    ci_lows[[k]] <- fit$ci_low
    ci_highs[[k]] <- fit$ci_high
  }
  out <- data.frame(
    label = labs,
    n_label = counts_cases + counts_controls,
    n_cases = counts_cases,
    n_controls = counts_controls,
    total_cases = total_cases,
    total_controls = total_controls,
    or = ors,
    ci_low = ci_lows,
    ci_high = ci_highs,
    p_value = pvals,
    stringsAsFactors = FALSE
  )
  out$q_value <- p.adjust(out$p_value, method = "BH")
  out[order(out$q_value, out$p_value, out$label), , drop = FALSE]
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

classify_bucket <- function(df) {
  cols <- intersect(c("phenotype_hint", "disease_subtype_hint", "final_case_control_label", "case_control_hint"), names(df))
  combined <- apply(df[, cols, drop = FALSE], 1, function(x) paste(tolower(as.character(x)), collapse = " "))
  out <- ifelse(
    grepl("healthy", combined), "healthy",
    ifelse(
      grepl("crohn|\\bccd\\b|\\bicd\\b", combined), "CD",
      ifelse(grepl("ulcerative|\\buc\\b", combined), "UC", "other")
    )
  )
  out
}

breakdown_string <- function(values) {
  sprintf("CD %d / UC %d / healthy %d", sum(values == "CD"), sum(values == "UC"), sum(values == "healthy"))
}

cohort_context_row <- function(cohort) {
  primary_comp <- switch(
    cohort,
    halfvarson_2017 = "IBD_vs_Healthy",
    hmp2 = "IBD_vs_Healthy",
    gevers_2014 = "Crohn_vs_Healthy",
    prjeb84421 = "OFG_vs_Healthy",
    jacobs_2023_ibs_250bp = "IBS_vs_Healthy"
  )
  primary <- load_assignment_csv(transfer_root, cohort, primary_comp)
  if (is.null(primary)) {
    return(NULL)
  }
  repeated_subjects <- sum(table(primary$subject_id) > 1L)
  subject_n <- length(unique(primary$subject_id))
  spec <- overview_specs[overview_specs$cohort == cohort, , drop = FALSE][1, ]

  if (cohort %in% c("halfvarson_2017", "hmp2")) {
    crohn <- load_assignment_csv(transfer_root, cohort, "Crohn_vs_Healthy")
    uc <- load_assignment_csv(transfer_root, cohort, "UC_vs_Healthy")
    crohn_cond <- derive_condition(crohn, "Crohn_vs_Healthy")
    uc_cond <- derive_condition(uc, "UC_vs_Healthy")
    sample_text <- sprintf(
      "%s %s samples (%d/%d); subset contrasts used %d/%d Crohn-vs-healthy and %d/%d UC-vs-healthy samples.",
      fmt_int(nrow(primary)),
      if (primary_comp == "IBD_vs_Healthy") "pooled IBD-vs-healthy" else comparison_display(primary_comp),
      sum(primary$condition <- derive_condition(primary, primary_comp) == 1L),
      sum(primary$condition == 0L),
      sum(crohn_cond == 1L), sum(crohn_cond == 0L),
      sum(uc_cond == 1L), sum(uc_cond == 0L)
    )
  } else {
    cond <- derive_condition(primary, primary_comp)
    sample_text <- sprintf(
      "%s %s samples (%d/%d).",
      fmt_int(nrow(primary)),
      if (primary_comp == "Crohn_vs_Healthy") "Crohn-vs-healthy" else if (primary_comp == "OFG_vs_Healthy") "OFG-vs-healthy" else if (primary_comp == "IBS_vs_Healthy") "IBS-vs-healthy baseline stool" else comparison_display(primary_comp),
      sum(cond == 1L),
      sum(cond == 0L)
    )
  }

  caveat <- if (repeated_subjects > 0L) {
    sprintf("%s; %s detected subjects, %s with repeated samples.", spec$context_prefix, fmt_int(subject_n), fmt_int(repeated_subjects))
  } else {
    sprintf("%s; no repeated subjects detected.", spec$context_prefix)
  }

  data.frame(
    Cohort = spec$cohort_short,
    samples_and_labels = sample_text,
    context_and_caveat = caveat,
    stringsAsFactors = FALSE
  )
}

render_table <- function(path, colspec, header, rows, size = "\\scriptsize") {
  lines <- c(
    "\\begin{center}",
    size,
    sprintf("\\begin{tabularx}{\\textwidth}{@{} %s @{}}", colspec),
    "\\toprule",
    paste(header, collapse = " & "),
    " \\\\",
    "\\midrule",
    rows,
    "\\bottomrule",
    "\\end{tabularx}",
    "\\end{center}"
  )
  writeLines(lines, path)
}

build_antibiotic_sensitivity <- function() {
  assign <- fread(cmd = sprintf("gzip -dc %s", shQuote(file.path(run_dir, "agp_absorb_assignments.tsv.gz"))), data.table = FALSE)
  meta <- fread(
    cmd = sprintf("gzip -dc %s", shQuote(agp_metadata_path)),
    select = c("Run", "BioProject", "Antibiotic_Status", "Antibiotics_Use"),
    data.table = FALSE
  )
  if (!"Antibiotic_Status" %in% names(meta) && "Antibiotics_Use" %in% names(meta)) {
    meta$Antibiotic_Status <- derive_agp_antibiotic_status(meta$Antibiotics_Use)
  }
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
        label = spec$label,
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
  fwrite(out, file.path(tables_dir, "TABLE_S8_ibd_antibiotic_status_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")

  analysis_order <- c("age_sex_bmi_complete_case", "plus_antibiotic_status", "excluding_week_month_6mo_antibiotics")
  render_rows <- character()
  for (i in seq_len(nrow(label_specs))) {
    spec <- label_specs[i, ]
    vals <- lapply(analysis_order, function(a) {
      row <- out[out$label == spec$label & out$analysis == a, , drop = FALSE]
      sprintf("OR %s; q = %s", fmt_or(row$adjusted_or[[1]]), fmt_prob(row$q_value_selected_labels[[1]]))
    })
    render_rows[[i]] <- sprintf(
      "%s & %d & %s & %s & %s \\\\",
      format_label_tex(spec$label),
      spec$depth,
      vals[[1]],
      vals[[2]],
      vals[[3]]
    )
  }
  render_table(
    file.path(tables_dir, "TABLE_S8_ibd_antibiotic_status_sensitivity.tex"),
    "Y r Y Y Y",
    c("Label", "Depth", "Age/sex/BMI", "Plus antibiotic status", "Excluding recent antibiotics"),
    render_rows
  )
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
        label = spec$label,
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
  fwrite(out, file.path(tables_dir, "TABLE_S9_absorb_policy_selected_ibd_robustness.tsv"), sep = "\t", quote = FALSE, na = "NA")

  render_rows <- character()
  for (i in seq_len(nrow(label_specs))) {
    spec <- label_specs[i, ]
    absorb <- out[out$label == spec$label & out$view == "absorb", , drop = FALSE]
    rare <- out[out$label == spec$label & out$view == "rare_policy", , drop = FALSE]
    render_rows[[i]] <- sprintf(
      "%s & %d & OR %s; q = %s & OR %s; q = %s \\\\",
      format_label_tex(spec$label),
      spec$depth,
      fmt_or(absorb$fisher_or[[1]]),
      fmt_prob(absorb$q_value_selected_labels[[1]]),
      fmt_or(rare$fisher_or[[1]]),
      fmt_prob(rare$q_value_selected_labels[[1]])
    )
  }
  render_table(
    file.path(tables_dir, "TABLE_S9_absorb_policy_selected_ibd_robustness.tex"),
    "Y r Y Y",
    c("Label", "Depth", "Absorb view", "Pure-policy baseline"),
    render_rows
  )
  out
}

build_taxonomic_rank_table <- function() {
  components <- unique(unlist(strsplit(label_specs$label, "__", fixed = TRUE)))
  rows <- lapply(components, function(component) {
    hit <- taxon_reference[taxon_reference$feature_label == component, , drop = FALSE]
    if (!nrow(hit)) {
      return(data.frame(component = component, taxonomic_rank = "not_found", taxonomy = "", stringsAsFactors = FALSE))
    }
    data.frame(component = component, taxonomic_rank = parse_rank(hit$taxonomy[[1]]), taxonomy = hit$taxonomy[[1]], stringsAsFactors = FALSE)
  })
  out <- rbindlist(rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S10_selected_ibd_label_taxonomic_ranks.tsv"), sep = "\t", quote = FALSE, na = "NA")

  render_rows <- vapply(seq_len(nrow(out)), function(i) {
    row <- out[i, ]
    sprintf(
      "%s & %s & %s \\\\",
      format_component_tex(row$component),
      tex_escape(paste0(toupper(substring(row$taxonomic_rank, 1L, 1L)), substring(row$taxonomic_rank, 2L))),
      tex_escape(summarize_taxonomy(row$taxonomy))
    )
  }, character(1L))
  render_table(
    file.path(tables_dir, "TABLE_S10_selected_ibd_label_taxonomic_ranks.tex"),
    "l l Y",
    c("Component", "Rank used by retained feature", "SILVA lineage summary"),
    render_rows,
    size = "\\footnotesize"
  )
  out
}

build_external_validation_overview <- function() {
  rebuilt_root <- file.path(validation_dir, "absorb")

  s11a_rows <- list()
  s11b_rows <- list()
  for (i in seq_len(nrow(overview_specs))) {
    spec <- overview_specs[i, ]
    rebuilt <- read_best_signal(rebuilt_root, spec$cohort, spec$comparison)
    transfer <- read_best_signal(transfer_root, spec$cohort, spec$comparison)
    mapping_path <- file.path(transfer_root, spec$cohort, sprintf("%s__%s_mapping_summary.csv", spec$cohort, spec$comparison))
    mapping <- if (file.exists(mapping_path)) read.csv(mapping_path, check.names = FALSE) else data.frame()
    depth2_cov <- if (nrow(mapping) && any(mapping$depth == 2L)) 100 * mapping$mapping_rate_total[mapping$depth == 2L][1] else NA_real_

    s11a_rows[[i]] <- data.frame(
      Cohort = spec$cohort_short,
      Comparison = comparison_display(spec$comparison),
      best_signal = if (!is.null(rebuilt)) sprintf("%s (depth %d, q = %s)", format_label_tex(rebuilt$DCST[[1]]), rebuilt$depth[[1]], fmt_prob(rebuilt$q_value[[1]])) else "NA",
      interpretation = interpret_overview("rebuilt", spec$comparison, rebuilt$q_value[[1]] %||% NA_real_),
      stringsAsFactors = FALSE
    )

    s11b_rows[[i]] <- data.frame(
      Cohort = spec$cohort_short,
      Comparison = comparison_display(spec$comparison),
      best_signal = if (!is.null(transfer)) sprintf("%s (depth %d, q = %s)", format_label_tex(transfer$DCST[[1]]), transfer$depth[[1]], fmt_prob(transfer$q_value[[1]])) else "NA",
      depth2_coverage = depth2_cov,
      interpretation = interpret_overview("transfer", spec$comparison, transfer$q_value[[1]] %||% NA_real_, depth2_cov),
      stringsAsFactors = FALSE
    )
  }
  s11a <- rbindlist(s11a_rows, fill = TRUE)
  s11b <- rbindlist(s11b_rows, fill = TRUE)
  s11c <- rbindlist(lapply(unique(overview_specs$cohort), cohort_context_row), fill = TRUE)

  s11d_rows <- list()
  for (i in seq_len(nrow(focus_specs))) {
    spec <- focus_specs[i, ]
    assign <- load_assignment_csv(transfer_root, spec$cohort, spec$comparison)
    if (is.null(assign)) next
    cond <- derive_condition(assign, spec$comparison)
    controls_sample <- sum(cond == 0L)
    controls_subject <- sum(derive_condition(assign[!duplicated(assign$subject_id), , drop = FALSE], spec$comparison) == 0L)
    sample_breakdown <- NULL
    subject_breakdown <- NULL
    if (spec$comparison == "Crohn_vs_Healthy") {
      sample_breakdown <- sprintf("CD %d / UC 0 / healthy %d", sum(cond == 1L), controls_sample)
      subject_breakdown <- sprintf("CD %d / UC 0 / healthy %d", sum(!duplicated(assign$subject_id) & cond == 1L), controls_subject)
    } else if (spec$comparison == "UC_vs_Healthy") {
      sample_breakdown <- sprintf("CD 0 / UC %d / healthy %d", sum(cond == 1L), controls_sample)
      subject_breakdown <- sprintf("CD 0 / UC %d / healthy %d", sum(!duplicated(assign$subject_id) & cond == 1L), controls_subject)
    } else if (spec$comparison == "IBD_vs_Healthy" && spec$cohort %in% c("halfvarson_2017", "hmp2")) {
      crohn <- load_assignment_csv(transfer_root, spec$cohort, "Crohn_vs_Healthy")
      uc <- load_assignment_csv(transfer_root, spec$cohort, "UC_vs_Healthy")
      crohn_cond <- derive_condition(crohn, "Crohn_vs_Healthy")
      uc_cond <- derive_condition(uc, "UC_vs_Healthy")
      sample_breakdown <- sprintf("CD %d / UC %d / healthy %d", sum(crohn_cond == 1L), sum(uc_cond == 1L), controls_sample)
      subject_breakdown <- sprintf(
        "CD %d / UC %d / healthy %d",
        sum(!duplicated(crohn$subject_id) & crohn_cond == 1L),
        sum(!duplicated(uc$subject_id) & uc_cond == 1L),
        controls_subject
      )
    } else {
      bucket_sample <- classify_bucket(assign)
      assign_subject <- assign[!duplicated(assign$subject_id), , drop = FALSE]
      bucket_subject <- classify_bucket(assign_subject)
      sample_breakdown <- breakdown_string(bucket_sample)
      subject_breakdown <- breakdown_string(bucket_subject)
    }
    repeated_subjects <- sum(table(assign$subject_id) > 1L)
    selection_strategy <- if (repeated_subjects > 0L) {
      "Primary sample-level analysis; first-sample and bootstrap subject-level sensitivity."
    } else {
      "All eligible filtered samples retained; repeated-subject burden is low."
    }
    if (spec$cohort == "gevers_2014") {
      selection_strategy <- "All eligible filtered samples retained; repeated-subject burden is low."
    }
    s11d_rows[[length(s11d_rows) + 1L]] <- data.frame(
      Cohort = spec$cohort_short,
      Comparison = comparison_display(spec$comparison),
      Samples = nrow(assign),
      Subjects = sprintf("%s (%s repeated)", fmt_int(length(unique(assign$subject_id))), fmt_int(repeated_subjects)),
      sample_breakdown = sample_breakdown,
      subject_breakdown = subject_breakdown,
      selection_strategy = selection_strategy,
      stringsAsFactors = FALSE
    )
  }
  s11d <- rbindlist(s11d_rows, fill = TRUE)

  fwrite(s11a, file.path(tables_dir, "TABLE_S11A_rebuilt_external_validation_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(s11b, file.path(tables_dir, "TABLE_S11B_transfer_external_validation_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(s11c, file.path(tables_dir, "TABLE_S11C_external_validation_cohort_context.tsv"), sep = "\t", quote = FALSE, na = "NA")
  fwrite(s11d, file.path(tables_dir, "TABLE_S11D_external_ibd_cohort_characterization.tsv"), sep = "\t", quote = FALSE, na = "NA")

  render_table(
    file.path(tables_dir, "TABLE_S11A_rebuilt_external_validation_summary.tex"),
    "l l Y Y",
    c("Cohort", "Comparison", "Best rebuilt signal", "Interpretation"),
    vapply(seq_len(nrow(s11a)), function(i) {
      row <- s11a[i, ]
      sprintf("%s & %s & %s & %s \\\\", tex_escape(row$Cohort), tex_escape(row$Comparison), row$best_signal, tex_escape(row$interpretation))
    }, character(1L)),
    size = "\\footnotesize"
  )

  render_table(
    file.path(tables_dir, "TABLE_S11B_transfer_external_validation_summary.tex"),
    "l l Y r Y",
    c("Cohort", "Comparison", "Best AGP-derived label-transfer signal", "Depth-2 coverage (\\%)", "Interpretation"),
    vapply(seq_len(nrow(s11b)), function(i) {
      row <- s11b[i, ]
      sprintf(
        "%s & %s & %s & %s & %s \\\\",
        tex_escape(row$Cohort),
        tex_escape(row$Comparison),
        row$best_signal,
        fmt_pct(row$depth2_coverage),
        tex_escape(row$interpretation)
      )
    }, character(1L)),
    size = "\\footnotesize"
  )

  render_table(
    file.path(tables_dir, "TABLE_S11C_external_validation_cohort_context.tex"),
    "l Y Y",
    c("Cohort", "Samples and comparison labels", "Context and caveat"),
    vapply(seq_len(nrow(s11c)), function(i) {
      row <- s11c[i, ]
      sprintf("%s & %s & %s \\\\", tex_escape(row$Cohort), tex_escape(row$samples_and_labels), tex_escape(row$context_and_caveat))
    }, character(1L))
  )

  render_table(
    file.path(tables_dir, "TABLE_S11D_external_ibd_cohort_characterization.tex"),
    "l l r r Y Y Y",
    c("Cohort", "Comparison", "Samples", "Subjects", "Sample breakdown (CD / UC / healthy)", "Subject breakdown (CD / UC / healthy)", "Selection strategy"),
    vapply(seq_len(nrow(s11d)), function(i) {
      row <- s11d[i, ]
      sprintf(
        "%s & %s & %s & %s & %s & %s & %s \\\\",
        tex_escape(row$Cohort),
        tex_escape(row$Comparison),
        fmt_int(row$Samples),
        tex_escape(row$Subjects),
        tex_escape(row$sample_breakdown),
        tex_escape(row$subject_breakdown),
        tex_escape(row$selection_strategy)
      )
    }, character(1L))
  )

  list(s11a = s11a, s11b = s11b, s11c = s11c, s11d = s11d)
}

evaluate_selected_label <- function(labels, condition, target_label) {
  mapped <- !is.na(labels) & nzchar(labels)
  labels_m <- as.character(labels[mapped])
  condition_m <- as.integer(condition[mapped])
  if (!length(labels_m) || length(unique(condition_m)) < 2L) {
    return(list(q_value = 1, corrected_or = NA_real_))
  }
  assoc <- association_table(labels_m, condition_m)
  target_present <- labels_m == target_label
  a <- sum(target_present & condition_m == 1L)
  b <- sum(target_present & condition_m == 0L)
  c <- sum(!target_present & condition_m == 1L)
  d <- sum(!target_present & condition_m == 0L)
  row <- assoc[assoc$label == target_label, , drop = FALSE]
  q_value <- if (nrow(row)) row$q_value[[1]] else 1
  list(q_value = q_value, corrected_or = haldane_or(a, b, c, d))
}

bootstrap_subject_consensus <- function(assign, label_col, target_label, condition, iterations, seed) {
  subject_groups <- split(seq_len(nrow(assign)), assign$subject_id)
  set.seed(seed)
  qs <- numeric(iterations)
  ors <- numeric(iterations)
  for (iter in seq_len(iterations)) {
    selected_idx <- vapply(subject_groups, function(idx) {
      if (length(idx) == 1L) idx[[1L]] else sample(idx, size = 1L)
    }, integer(1L))
    selected <- assign[selected_idx, , drop = FALSE]
    eval <- evaluate_selected_label(selected[[label_col]], condition[selected_idx], target_label)
    qs[[iter]] <- eval$q_value
    ors[[iter]] <- eval$corrected_or
  }
  finite_or <- ors[is.finite(ors)]
  list(
    iterations = iterations,
    consensus_rate = mean(qs < 0.05, na.rm = TRUE),
    median_q = median(qs, na.rm = TRUE),
    median_or = if (length(finite_or)) median(finite_or, na.rm = TRUE) else NA_real_,
    or_lo = if (length(finite_or)) as.numeric(quantile(finite_or, probs = 0.025, na.rm = TRUE)) else NA_real_,
    or_hi = if (length(finite_or)) as.numeric(quantile(finite_or, probs = 0.975, na.rm = TRUE)) else NA_real_
  )
}

order_samples_within_subject <- function(df) {
  timepoint <- if ("timepoint_hint" %in% names(df)) {
    suppressWarnings(as.numeric(df$timepoint_hint))
  } else {
    rep(NA_real_, nrow(df))
  }
  ord <- order(df$subject_id, timepoint, df$run_accession, na.last = TRUE)
  df[ord, , drop = FALSE]
}

first_subject_subset <- function(df) {
  ordered <- order_samples_within_subject(df)
  ordered[!duplicated(ordered$subject_id), , drop = FALSE]
}

build_subject_level_diagnostics <- function() {
  modes <- data.frame(
    mode = c("rebuilt_cohort", "agp_label_transfer"),
    root = c(file.path(validation_dir, "absorb"), transfer_root),
    label_template = c("depth%d_absorb", "depth%d_frozen_agp_absorb"),
    stringsAsFactors = FALSE
  )

  rows <- list()
  row_id <- 0L
  for (i in seq_len(nrow(subject_diag_specs))) {
    spec <- subject_diag_specs[i, ]
    combined <- list()
    for (j in seq_len(nrow(modes))) {
      row_id <- row_id + 1L
      mode <- modes[j, ]
      best <- read_best_signal(mode$root, spec$cohort, spec$comparison)
      if (is.null(best)) next
      assign <- load_assignment_csv(mode$root, spec$cohort, spec$comparison)
      if (is.null(assign)) next
      assign$condition <- derive_condition(assign, spec$comparison)
      assign <- assign[order(assign$subject_id, assign$run_accession), , drop = FALSE]
      first <- assign[!duplicated(assign$subject_id), , drop = FALSE]
      label_col <- sprintf(mode$label_template, best$depth[[1]])
      first_eval <- evaluate_selected_label(first[[label_col]], first$condition, as.character(best$DCST[[1]]))
      boot <- bootstrap_subject_consensus(
        assign = assign,
        label_col = label_col,
        target_label = as.character(best$DCST[[1]]),
        condition = assign$condition,
        iterations = bootstrap_iterations,
        seed = bootstrap_seed + row_id
      )
      combined[[mode$mode]] <- data.frame(
        cohort = spec$cohort_short,
        comparison = comparison_display(spec$comparison),
        mode = mode$mode,
        depth = best$depth[[1]],
        sample_best_label = as.character(best$DCST[[1]]),
        sample_q = best$q_value[[1]],
        first_subject_q = first_eval$q_value,
        bootstrap_iterations = boot$iterations,
        bootstrap_consensus_rate = boot$consensus_rate,
        bootstrap_median_q = boot$median_q,
        bootstrap_median_or = boot$median_or,
        bootstrap_or_lo = boot$or_lo,
        bootstrap_or_hi = boot$or_hi,
        n_samples = nrow(assign),
        n_subjects = length(unique(assign$subject_id)),
        repeated_subjects = sum(table(assign$subject_id) > 1L),
        stringsAsFactors = FALSE
      )
    }
    if (length(combined)) {
      rows[[length(rows) + 1L]] <- do.call(rbind, combined)
    }
  }
  out <- rbindlist(rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S12_external_subject_level_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")

  as_cell <- function(row) {
    sprintf(
      "%s\\;(<0.05); med q %s; OR %s [%s, %s]",
      fmt_pct(100 * row$bootstrap_consensus_rate),
      fmt_prob(row$bootstrap_median_q),
      fmt_or(row$bootstrap_median_or),
      fmt_or(row$bootstrap_or_lo),
      fmt_or(row$bootstrap_or_hi)
    )
  }

  display_rows <- character()
  unique_pairs <- unique(out[, c("cohort", "comparison")])
  for (i in seq_len(nrow(unique_pairs))) {
    pair <- unique_pairs[i, ]
    rebuilt <- out[out$cohort == pair$cohort & out$comparison == pair$comparison & out$mode == "rebuilt_cohort", , drop = FALSE]
    transfer <- out[out$cohort == pair$cohort & out$comparison == pair$comparison & out$mode == "agp_label_transfer", , drop = FALSE]
    if (!nrow(rebuilt) || !nrow(transfer)) next
    display_rows[[length(display_rows) + 1L]] <- sprintf(
      "%s & %s & %s $\\rightarrow$ %s & %s $\\rightarrow$ %s & %s & %s \\\\",
      tex_escape(pair$cohort),
      tex_escape(pair$comparison),
      fmt_prob(rebuilt$sample_q[[1]]),
      fmt_prob(rebuilt$first_subject_q[[1]]),
      fmt_prob(transfer$sample_q[[1]]),
      fmt_prob(transfer$first_subject_q[[1]]),
      as_cell(rebuilt),
      as_cell(transfer)
    )
  }

  render_table(
    file.path(tables_dir, "TABLE_S12_external_subject_level_sensitivity.tex"),
    "l l Y Y Y Y",
    c("Cohort", "Comparison", "Rebuilt q (sample $\\rightarrow$ first subject)", "AGP-transfer q (sample $\\rightarrow$ first subject)", "Rebuilt bootstrap consensus", "AGP-transfer bootstrap consensus"),
    display_rows
  )
  out
}

load_halfvarson_clinical_metadata <- function() {
  validation_meta <- fread(
    file.path(gut_root, "outputs", "validation_ingest", "external_harmonized", "halfvarson_2017", "validation_metadata_table.tsv"),
    data.table = FALSE
  )
  clinical_meta <- fread(
    file.path(gut_root, "data", "validation_cohorts", "external", "halfvarson_2017", "metadata", "halfvarson_supplementary_dataset1_parsed.tsv"),
    data.table = FALSE
  )
  validation_meta$participant_id_hint <- as.character(validation_meta$participant_id_hint)
  validation_meta$timepoint_hint <- suppressWarnings(as.numeric(validation_meta$timepoint_hint))
  clinical_meta$host_subject_id <- as.character(clinical_meta$host_subject_id)
  clinical_meta$timepoint <- suppressWarnings(as.numeric(clinical_meta$timepoint))
  merge(
    validation_meta,
    clinical_meta,
    by.x = c("participant_id_hint", "timepoint_hint"),
    by.y = c("host_subject_id", "timepoint"),
    all.x = TRUE
  )
}

load_gevers_medication_metadata <- function() {
  validation_meta <- fread(
    file.path(gut_root, "outputs", "validation_ingest", "external_harmonized", "gevers_2014", "gevers_2014_validation_metadata.tsv"),
    data.table = FALSE
  )
  attr_paths <- c(
    file.path(gut_root, "data", "validation_cohorts", "external", "gevers_2014", "metadata", "ena_sample_attributes_PRJEB13679.tsv"),
    file.path(gut_root, "data", "validation_cohorts", "external", "gevers_2014", "metadata", "ena_sample_attributes_PRJEB13680.tsv")
  )
  attrs <- rbindlist(lapply(attr_paths, function(path) fread(path, data.table = FALSE)), fill = TRUE)
  attrs <- attrs[!duplicated(attrs$sample_accession), , drop = FALSE]
  merged <- merge(
    validation_meta,
    attrs[, c("sample_accession", "antibiotics", "mesalamine", "steroids", "biologics"), drop = FALSE],
    by = "sample_accession",
    all.x = TRUE
  )
  for (field in c("antibiotics", "mesalamine", "steroids", "biologics")) {
    merged[[field]] <- normalize_binary_flag(merged[[field]])
  }
  merged[, c("run_accession", "sample_accession", "participant_id_hint", "analysis_label", "antibiotics", "mesalamine", "steroids", "biologics")]
}

build_halfvarson_clinical_followup <- function() {
  clinical_meta <- load_halfvarson_clinical_metadata()
  cd_codes <- c("CCD", "CC", "ICD_r", "ICD_nr", "LC")

  analysis_specs <- list(
    list(
      table_id = "S13A",
      key = "uc_extent",
      title = "Halfvarson UC extent follow-up",
      comparison = "UC_vs_Healthy",
      summary_sentence = "Cases were extensive UC samples (E3) and controls were distal UC samples (E1/E2).",
      subset_fn = function(df) df$ibd_subtype == "UC" & df$uc_extent %in% c("Left sided (E2)", "Proctitis (E1)", "Extensive (E3)"),
      condition_fn = function(df) as.integer(df$uc_extent == "Extensive (E3)"),
      positive_label = "Extensive UC (E3)",
      negative_label = "Distal UC (E1/E2)",
      context_fn = function(df) {
        group <- ifelse(df$uc_extent == "Extensive (E3)", "Extensive UC (E3)", "Distal UC (E1/E2)")
        subject_group <- tapply(group, df$subject_id, function(x) unique(x)[1L])
        sprintf(
          "The sample-level follow-up included %s extensive and %s distal samples from %s and %s subjects, respectively.",
          fmt_int(sum(group == "Extensive UC (E3)")),
          fmt_int(sum(group == "Distal UC (E1/E2)")),
          fmt_int(sum(subject_group == "Extensive UC (E3)")),
          fmt_int(sum(subject_group == "Distal UC (E1/E2)"))
        )
      }
    ),
    list(
      table_id = "S13B",
      key = "cd_location",
      title = "Halfvarson Crohn disease location follow-up",
      comparison = "Crohn_vs_Healthy",
      summary_sentence = "Cases were ileal-involved Crohn samples (Montreal L1/L1+L4 and L3/L3+L4) and controls were colonic-only Crohn samples (Montreal L2).",
      subset_fn = function(df) {
        df$ibd_subtype %in% cd_codes &
          df$cd_location %in% c("Colonic (L2)", "Ileal (L1)", "Ileocolonic (L3)", "Ileal and Upper-GI (L1+L4)", "Ileocolonic and Upper-GI (L3+L4)")
      },
      condition_fn = function(df) as.integer(df$cd_location %in% c("Ileal (L1)", "Ileocolonic (L3)", "Ileal and Upper-GI (L1+L4)", "Ileocolonic and Upper-GI (L3+L4)")),
      positive_label = "Ileal-involved Crohn disease",
      negative_label = "Colonic-only Crohn disease",
      context_fn = function(df) {
        group <- ifelse(
          df$cd_location %in% c("Ileal (L1)", "Ileocolonic (L3)", "Ileal and Upper-GI (L1+L4)", "Ileocolonic and Upper-GI (L3+L4)"),
          "Ileal-involved Crohn disease",
          "Colonic-only Crohn disease"
        )
        subject_group <- tapply(group, df$subject_id, function(x) unique(x)[1L])
        sprintf(
          "The sample-level follow-up included %s ileal-involved and %s colonic-only samples from %s and %s subjects, respectively.",
          fmt_int(sum(group == "Ileal-involved Crohn disease")),
          fmt_int(sum(group == "Colonic-only Crohn disease")),
          fmt_int(sum(subject_group == "Ileal-involved Crohn disease")),
          fmt_int(sum(subject_group == "Colonic-only Crohn disease"))
        )
      }
    ),
    list(
      table_id = "S13C",
      key = "cd_calprotectin",
      title = "Halfvarson Crohn disease calprotectin follow-up",
      comparison = "Crohn_vs_Healthy",
      summary_sentence = "Cases were Crohn samples with fecal calprotectin $\\geq 250$ \\textmu g/g and controls were Crohn samples with fecal calprotectin $< 250$ \\textmu g/g; this was used as a conservative high-inflammatory-burden split rather than a remission cut-off.",
      subset_fn = function(df) df$ibd_subtype %in% cd_codes & is.finite(df$cal_num),
      condition_fn = function(df) as.integer(df$cal_num >= 250),
      positive_label = "Crohn calprotectin $\\geq 250$ \\textmu g/g",
      negative_label = "Crohn calprotectin $< 250$ \\textmu g/g",
      context_fn = function(df) {
        group <- ifelse(df$cal_num >= 250, "High", "Low")
        by_subject <- tapply(group, df$subject_id, function(x) paste(sort(unique(x)), collapse = "+"))
        sprintf(
          "The sample-level follow-up included %s high- and %s low-calprotectin samples across %s Crohn subjects; %s subjects contributed both high- and low-calprotectin samples over time.",
          fmt_int(sum(group == "High")),
          fmt_int(sum(group == "Low")),
          fmt_int(length(unique(df$subject_id))),
          fmt_int(sum(by_subject == "High+Low"))
        )
      }
    )
  )

  mode_specs <- list(
    list(mode = "rebuilt", display = "Rebuilt", root = file.path(validation_dir, "absorb"), label_template = "depth%d_absorb"),
    list(mode = "agp_transfer", display = "AGP transfer", root = transfer_root, label_template = "depth%d_frozen_agp_absorb")
  )

  tex_lines <- character()
  summary_rows <- list()
  subject_sensitivity_rows <- list()
  subject_row_id <- 0L

  for (spec in analysis_specs) {
    raw_rows <- list()
    display_rows <- character()
    context_text <- NULL

    for (mode_spec in mode_specs) {
      assign <- load_assignment_csv(mode_spec$root, "halfvarson_2017", spec$comparison)
      if (is.null(assign)) {
        next
      }
      assign <- merge(
        assign,
        clinical_meta[, c("run_accession", "timepoint_hint", "ibd_subtype", "uc_extent", "cd_location", "calprotectin")],
        by = "run_accession",
        all.x = TRUE
      )
      assign$cal_num <- suppressWarnings(as.numeric(assign$calprotectin))
      keep <- spec$subset_fn(assign)
      assign <- assign[keep, , drop = FALSE]
      assign <- order_samples_within_subject(assign)
      condition <- spec$condition_fn(assign)
      if (is.null(context_text)) {
        context_text <- spec$context_fn(assign)
      }

      for (depth in 1:4) {
        label_col <- sprintf(mode_spec$label_template, depth)
        assoc <- association_table(assign[[label_col]], condition)
        if (!nrow(assoc)) {
          next
        }
        assoc$mode <- mode_spec$display
        assoc$depth <- depth
        assoc$analysis <- spec$key
        raw_rows[[length(raw_rows) + 1L]] <- assoc

        top <- assoc[1, , drop = FALSE]
        first <- first_subject_subset(assign)
        first_condition <- spec$condition_fn(first)
        first_eval <- evaluate_selected_label(first[[label_col]], first_condition, top$label[[1]])
        subject_row_id <- subject_row_id + 1L
        boot <- bootstrap_subject_consensus(
          assign = assign,
          label_col = label_col,
          target_label = top$label[[1]],
          condition = condition,
          iterations = bootstrap_iterations,
          seed = bootstrap_seed + 10000L + subject_row_id
        )
        a <- top$n_cases[[1]]
        b <- top$n_controls[[1]]
        c <- top$total_cases[[1]] - a
        d <- top$total_controls[[1]] - b
        bayes <- bayes_or_ci(a, b, c, d, sprintf("%s|%s|%s|%d|%s", spec$key, mode_spec$mode, spec$comparison, depth, top$label[[1]]))
        summary_rows[[length(summary_rows) + 1L]] <- data.frame(
          analysis = spec$key,
          table_id = spec$table_id,
          mode = mode_spec$display,
          depth = depth,
          positive_samples = top$total_cases[[1]],
          negative_samples = top$total_controls[[1]],
          top_label = top$label[[1]],
          n_cases = a,
          n_controls = b,
          or = top$or[[1]],
          ci_low = top$ci_low[[1]],
          ci_high = top$ci_high[[1]],
          p_value = top$p_value[[1]],
          q_value = top$q_value[[1]],
          bayes_or = bayes$or,
          bayes_lo = bayes$lo,
          bayes_hi = bayes$hi,
          stringsAsFactors = FALSE
        )
        subject_sensitivity_rows[[length(subject_sensitivity_rows) + 1L]] <- data.frame(
          analysis = spec$key,
          table_id = spec$table_id,
          mode = mode_spec$display,
          depth = depth,
          top_label = top$label[[1]],
          sample_q = top$q_value[[1]],
          first_subject_q = first_eval$q_value,
          first_subject_positive = sum(first_condition == 1L, na.rm = TRUE),
          first_subject_negative = sum(first_condition == 0L, na.rm = TRUE),
          bootstrap_iterations = boot$iterations,
          bootstrap_consensus_rate = boot$consensus_rate,
          bootstrap_median_q = boot$median_q,
          bootstrap_median_or = boot$median_or,
          bootstrap_or_lo = boot$or_lo,
          bootstrap_or_hi = boot$or_hi,
          unique_subjects = length(unique(assign$subject_id)),
          repeated_subjects = sum(table(assign$subject_id) > 1L),
          stringsAsFactors = FALSE
        )
        display_rows[[length(display_rows) + 1L]] <- sprintf(
          "%s & %d & %s/%s & %s & OR %s [%s, %s] & %s & %s & OR %s [%s, %s] \\\\",
          tex_escape(mode_spec$display),
          depth,
          fmt_int(top$total_cases[[1]]),
          fmt_int(top$total_controls[[1]]),
          format_label_tex(top$label[[1]]),
          fmt_or(top$or[[1]]),
          fmt_or(top$ci_low[[1]]),
          fmt_or(top$ci_high[[1]]),
          fmt_prob(top$p_value[[1]]),
          fmt_prob(top$q_value[[1]]),
          fmt_or(bayes$or),
          fmt_or(bayes$lo),
          fmt_or(bayes$hi)
        )
      }
    }

    raw <- rbindlist(raw_rows, fill = TRUE)
    fwrite(raw, file.path(tables_dir, sprintf("TABLE_%s_halfvarson_%s.tsv", spec$table_id, spec$key)), sep = "\t", quote = FALSE, na = "NA")

    tex_lines <- c(
      tex_lines,
      sprintf("\\noindent\\textbf{Supplementary Table %s.} %s. %s %s Each row shows the top q-ranked dCST at the specified depth in rebuilt-cohort or direct AGP-derived transfer mode.", spec$table_id, spec$title, spec$summary_sentence, context_text),
      "",
      "\\begin{center}",
      "\\scriptsize",
      "\\begin{tabularx}{\\textwidth}{@{} l r r Y Y r r Y @{}}",
      "\\toprule",
      "Mode & Depth & Positive/negative samples & Top dCST & Frequentist OR (95\\% CI) & \\(p\\) & \\(q\\) & Bayesian OR (95\\% CrI) \\\\",
      "\\midrule",
      display_rows,
      "\\bottomrule",
      "\\end{tabularx}",
      "\\end{center}",
      "\\bigskip",
      ""
    )
  }

  summary <- rbindlist(summary_rows, fill = TRUE)
  fwrite(summary, file.path(tables_dir, "TABLE_S13A_S13C_halfvarson_clinical_followup_summary.tsv"), sep = "\t", quote = FALSE, na = "NA")
  writeLines(tex_lines, file.path(tables_dir, "TABLE_S13A_S13C_halfvarson_clinical_followup.tex"))
  subject_summary <- rbindlist(subject_sensitivity_rows, fill = TRUE)
  fwrite(subject_summary, file.path(tables_dir, "TABLE_S13D_halfvarson_clinical_followup_subject_sensitivity.tsv"), sep = "\t", quote = FALSE, na = "NA")

  analysis_display <- c(
    uc_extent = "UC extent",
    cd_location = "Crohn location",
    cd_calprotectin = "Crohn calprotectin"
  )

  subject_rows_tex <- vapply(seq_len(nrow(subject_summary)), function(i) {
    row <- subject_summary[i, ]
    bootstrap_cell <- sprintf(
      "%s\\;(<0.05); med q %s; OR %s [%s, %s]",
      fmt_pct(100 * row$bootstrap_consensus_rate),
      fmt_prob(row$bootstrap_median_q),
      fmt_or(row$bootstrap_median_or),
      fmt_or(row$bootstrap_or_lo),
      fmt_or(row$bootstrap_or_hi)
    )
    sprintf(
      "%s & %s & %d & %s & %s $\\rightarrow$ %s & %s \\\\",
      tex_escape(analysis_display[[row$analysis]] %||% row$analysis),
      tex_escape(row$mode),
      row$depth,
      format_label_tex(row$top_label),
      fmt_prob(row$sample_q),
      fmt_prob(row$first_subject_q),
      bootstrap_cell
    )
  }, character(1L))

  subject_tex <- c(
    "\\noindent\\textbf{Supplementary Table S13D.} Subject-aware sensitivity for the Halfvarson clinical stratification follow-up. For UC extent and Crohn location, the case/control definition was fixed at the subject level and the qualifier columns summarize one-sample-per-subject sensitivity for the sample-level top label from Supplementary Tables S13A and S13B. For Crohn calprotectin, the primary analysis remained visit-level; the qualifier columns therefore summarize one-sample-per-subject sensitivity for the same sample-level top label, with high-versus-low status re-evaluated from the selected sample in each resample. The bootstrap column reports the percentage of one-sample-per-subject resamples with $q < 0.05$, together with the median adjusted q-value and the median Jeffreys-corrected odds ratio with its empirical 95\\% percentile interval across \\(B = 1000\\) resamples.",
    "",
    "\\scriptsize",
    "\\setlength{\\LTleft}{0pt}",
    "\\setlength{\\LTright}{0pt}",
    "\\begin{longtable}{@{} p{0.12\\textwidth} p{0.11\\textwidth} p{0.05\\textwidth} p{0.25\\textwidth} p{0.17\\textwidth} p{0.24\\textwidth} @{}}",
    "\\toprule",
    "Analysis & Mode & Depth & Sample-level top dCST & Sample $\\rightarrow$ first-subject \\(q\\) & Bootstrap one-sample-per-subject consensus \\\\",
    "\\midrule",
    "\\endfirsthead",
    "\\toprule",
    "Analysis & Mode & Depth & Sample-level top dCST & Sample $\\rightarrow$ first-subject \\(q\\) & Bootstrap one-sample-per-subject consensus \\\\",
    "\\midrule",
    "\\endhead",
    "\\midrule",
    "\\multicolumn{6}{r}{\\emph{Continued on next page}} \\\\",
    "\\midrule",
    "\\endfoot",
    "\\bottomrule",
    "\\endlastfoot",
    subject_rows_tex,
    "\\end{longtable}",
    "\\normalsize"
  )
  writeLines(subject_tex, file.path(tables_dir, "TABLE_S13D_halfvarson_clinical_followup_subject_sensitivity.tex"))

  list(summary = summary, subject_summary = subject_summary)
}

build_gevers_medication_followup <- function() {
  gevers_meta <- load_gevers_medication_metadata()

  analysis_specs <- list(
    list(
      key = "abx_negative_cd_vs_healthy",
      medication_field = "antibiotics",
      title = "Antibiotic-negative Gevers subset",
      summary_sentence = "The cohort was restricted to samples explicitly flagged antibiotic-negative, and Crohn disease was compared with healthy controls within that subset.",
      subset_fn = function(df) df$antibiotics == "false",
      condition_fn = function(df) derive_condition(df, "Crohn_vs_Healthy"),
      context_fn = function(df, condition) sprintf(
        "The explicit antibiotic-negative subset contained %s Crohn samples and %s healthy-control samples.",
        fmt_int(sum(condition == 1L)),
        fmt_int(sum(condition == 0L))
      )
    ),
    list(
      key = "mesalamine_negative_cd_vs_healthy",
      medication_field = "mesalamine",
      title = "Mesalamine-negative Gevers subset",
      summary_sentence = "The cohort was restricted to samples explicitly flagged mesalamine-negative, and Crohn disease was compared with healthy controls within that subset.",
      subset_fn = function(df) df$mesalamine == "false",
      condition_fn = function(df) derive_condition(df, "Crohn_vs_Healthy"),
      context_fn = function(df, condition) sprintf(
        "The explicit mesalamine-negative subset contained %s Crohn samples and %s healthy-control samples.",
        fmt_int(sum(condition == 1L)),
        fmt_int(sum(condition == 0L))
      )
    ),
    list(
      key = "steroid_negative_cd_vs_healthy",
      medication_field = "steroids",
      title = "Steroid-negative Gevers subset",
      summary_sentence = "The cohort was restricted to samples explicitly flagged steroid-negative, and Crohn disease was compared with healthy controls within that subset.",
      subset_fn = function(df) df$steroids == "false",
      condition_fn = function(df) derive_condition(df, "Crohn_vs_Healthy"),
      context_fn = function(df, condition) sprintf(
        "The explicit steroid-negative subset contained %s Crohn samples and %s healthy-control samples.",
        fmt_int(sum(condition == 1L)),
        fmt_int(sum(condition == 0L))
      )
    ),
    list(
      key = "abx_positive_vs_negative_cd",
      medication_field = "antibiotics",
      title = "Antibiotic-positive vs antibiotic-negative Crohn samples",
      summary_sentence = "Only Crohn-case samples with known antibiotic status were retained, and antibiotic-positive samples were compared with antibiotic-negative samples.",
      subset_fn = function(df) derive_condition(df, "Crohn_vs_Healthy") == 1L & df$antibiotics %in% c("true", "false"),
      condition_fn = function(df) as.integer(df$antibiotics == "true"),
      context_fn = function(df, condition) sprintf(
        "Among Crohn-case samples with known antibiotic status, %s were antibiotic-positive and %s were antibiotic-negative.",
        fmt_int(sum(condition == 1L)),
        fmt_int(sum(condition == 0L))
      )
    ),
    list(
      key = "mesalamine_positive_vs_negative_cd",
      medication_field = "mesalamine",
      title = "Mesalamine-positive vs mesalamine-negative Crohn samples",
      summary_sentence = "Only Crohn-case samples with known mesalamine status were retained, and mesalamine-positive samples were compared with mesalamine-negative samples.",
      subset_fn = function(df) derive_condition(df, "Crohn_vs_Healthy") == 1L & df$mesalamine %in% c("true", "false"),
      condition_fn = function(df) as.integer(df$mesalamine == "true"),
      context_fn = function(df, condition) sprintf(
        "Among Crohn-case samples with known mesalamine status, %s were mesalamine-positive and %s were mesalamine-negative.",
        fmt_int(sum(condition == 1L)),
        fmt_int(sum(condition == 0L))
      )
    ),
    list(
      key = "steroid_positive_vs_negative_cd",
      medication_field = "steroids",
      title = "Steroid-positive vs steroid-negative Crohn samples",
      summary_sentence = "Only Crohn-case samples with known steroid status were retained, and steroid-positive samples were compared with steroid-negative samples.",
      subset_fn = function(df) derive_condition(df, "Crohn_vs_Healthy") == 1L & df$steroids %in% c("true", "false"),
      condition_fn = function(df) as.integer(df$steroids == "true"),
      context_fn = function(df, condition) sprintf(
        "Among Crohn-case samples with known steroid status, %s were steroid-positive and %s were steroid-negative.",
        fmt_int(sum(condition == 1L)),
        fmt_int(sum(condition == 0L))
      )
    )
  )

  mode_specs <- list(
    list(mode = "rebuilt", display = "Rebuilt", root = file.path(validation_dir, "absorb"), label_template = "depth%d_absorb"),
    list(mode = "agp_transfer", display = "AGP transfer", root = transfer_root, label_template = "depth%d_frozen_agp_absorb")
  )

  raw_rows <- list()
  tex_rows <- character()

  for (spec in analysis_specs) {
    context_text <- NULL
    for (mode_spec in mode_specs) {
      assign <- load_assignment_csv(mode_spec$root, "gevers_2014", "Crohn_vs_Healthy")
      if (is.null(assign)) {
        next
      }
      assign <- merge(
        assign,
        gevers_meta[, c("run_accession", "antibiotics", "mesalamine", "steroids", "biologics"), drop = FALSE],
        by = "run_accession",
        all.x = TRUE
      )
      keep <- spec$subset_fn(assign)
      assign <- assign[keep, , drop = FALSE]
      condition <- spec$condition_fn(assign)
      if (!nrow(assign) || length(unique(condition)) < 2L) {
        next
      }
      if (is.null(context_text)) {
        context_text <- spec$context_fn(assign, condition)
        tex_rows <- c(
          tex_rows,
          sprintf("\\noindent\\textbf{Analysis:} %s. %s %s", tex_escape(spec$title), spec$summary_sentence, context_text),
          "",
          "\\begin{center}",
          "\\scriptsize",
          "\\begin{tabularx}{\\textwidth}{@{} l r r Y Y r r Y @{}}",
          "\\toprule",
          "Mode & Depth & Positive/negative samples & Top dCST & Frequentist OR (95\\% CI) & \\(p\\) & \\(q\\) & Bayesian OR (95\\% CrI) \\\\",
          "\\midrule"
        )
      }

      for (depth in 1:4) {
        label_col <- sprintf(mode_spec$label_template, depth)
        assoc <- association_table(assign[[label_col]], condition)
        if (!nrow(assoc)) {
          next
        }
        top <- assoc[1, , drop = FALSE]
        a <- top$n_cases[[1]]
        b <- top$n_controls[[1]]
        c <- top$total_cases[[1]] - a
        d <- top$total_controls[[1]] - b
        bayes <- bayes_or_ci(
          a, b, c, d,
          sprintf("gevers_med|%s|%s|%s|%d|%s", spec$key, mode_spec$mode, spec$medication_field, depth, top$label[[1]])
        )

        raw_rows[[length(raw_rows) + 1L]] <- data.frame(
          analysis = spec$key,
          analysis_title = spec$title,
          medication_field = spec$medication_field,
          mode = mode_spec$display,
          depth = depth,
          positive_samples = top$total_cases[[1]],
          negative_samples = top$total_controls[[1]],
          top_label = top$label[[1]],
          n_positive_in_label = a,
          n_negative_in_label = b,
          or = top$or[[1]],
          ci_low = top$ci_low[[1]],
          ci_high = top$ci_high[[1]],
          p_value = top$p_value[[1]],
          q_value = top$q_value[[1]],
          bayes_or = bayes$or,
          bayes_lo = bayes$lo,
          bayes_hi = bayes$hi,
          stringsAsFactors = FALSE
        )

        tex_rows[[length(tex_rows) + 1L]] <- sprintf(
          "%s & %d & %s/%s & %s & OR %s [%s, %s] & %s & %s & OR %s [%s, %s] \\\\",
          tex_escape(mode_spec$display),
          depth,
          fmt_int(top$total_cases[[1]]),
          fmt_int(top$total_controls[[1]]),
          format_label_tex(top$label[[1]]),
          fmt_or(top$or[[1]]),
          fmt_or(top$ci_low[[1]]),
          fmt_or(top$ci_high[[1]]),
          fmt_prob(top$p_value[[1]]),
          fmt_prob(top$q_value[[1]]),
          fmt_or(bayes$or),
          fmt_or(bayes$lo),
          fmt_or(bayes$hi)
        )
      }
    }
    tex_rows <- c(
      tex_rows,
      "\\bottomrule",
      "\\end{tabularx}",
      "\\end{center}",
      "\\bigskip",
      ""
    )
  }

  out <- rbindlist(raw_rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S13E_gevers_medication_followup.tsv"), sep = "\t", quote = FALSE, na = "NA")

  tex_lines <- c(
    "\\noindent\\textbf{Supplementary Table S13E.} Exploratory Gevers medication follow-up. The first three analyses asked whether the pediatric Crohn portability boundary persisted after restricting the Gevers stool branch to samples explicitly flagged antibiotic-negative, mesalamine-negative, or steroid-negative. The next three analyses compared medication-positive with medication-negative Crohn-case samples among cases with known status for the same field. Biologics were not evaluable because no analyzed Gevers stool samples were flagged biologics-positive.",
    "",
    tex_rows
  )
  writeLines(tex_lines, file.path(tables_dir, "TABLE_S13E_gevers_medication_followup.tex"))
  out
}

build_mapping_balance <- function() {
  rows <- list()
  for (i in seq_len(nrow(mapping_diag_specs))) {
    spec <- mapping_diag_specs[i, ]
    path <- file.path(transfer_root, spec$cohort, sprintf("%s__%s_mapping_summary.csv", spec$cohort, spec$comparison))
    if (!file.exists(path)) next
    df <- read.csv(path, check.names = FALSE)
    row <- df[df$depth == spec$depth, , drop = FALSE]
    if (!nrow(row)) next
    rows[[length(rows) + 1L]] <- data.frame(
      cohort = spec$cohort_short,
      comparison = comparison_display(spec$comparison),
      depth = spec$depth,
      mapped_total = row$mapped_total[[1]],
      mapped_cases = row$mapped_cases[[1]],
      mapped_controls = row$mapped_controls[[1]],
      mapping_rate_total = row$mapping_rate_total[[1]],
      mapping_rate_cases = row$mapping_rate_cases[[1]],
      mapping_rate_controls = row$mapping_rate_controls[[1]],
      mapped_vs_unmapped_p = row$mapped_vs_unmapped_p[[1]],
      interpretation = interpret_mapping_balance(100 * row$mapping_rate_cases[[1]], 100 * row$mapping_rate_controls[[1]], row$mapped_vs_unmapped_p[[1]], spec$depth),
      stringsAsFactors = FALSE
    )
  }
  out <- rbindlist(rows, fill = TRUE)
  fwrite(out, file.path(tables_dir, "TABLE_S13_external_label_transfer_mapping_balance.tsv"), sep = "\t", quote = FALSE, na = "NA")

  render_table(
    file.path(tables_dir, "TABLE_S13_external_label_transfer_mapping_balance.tex"),
    "l r r r r Y",
    c(
      "Cohort comparison", "Depth", "Case coverage (\\%)", "Control coverage (\\%)", "Fisher p", "Note"
    ),
    vapply(seq_len(nrow(out)), function(i) {
      row <- out[i, ]
      sprintf(
        "%s %s & %d & %s & %s & %s & %s \\\\",
        tex_escape(row$cohort),
        tex_escape(row$comparison),
        row$depth,
        fmt_pct(100 * row$mapping_rate_cases),
        fmt_pct(100 * row$mapping_rate_controls),
        fmt_prob(row$mapped_vs_unmapped_p),
        tex_escape(row$interpretation)
      )
    }, character(1L)),
    size = "\\tiny"
  )
  out
}

write_summary <- function(cleaner_control, profiles, antibiotic, overview, subject, mapping, policy, ranks, clinical_followup, gevers_medication) {
  lines <- c(
    "# IBD Reviewer Sensitivity Asset Summary",
    "",
    "Generated by `build_ibd_reviewer_sensitivity_assets.R`.",
    "",
    sprintf("Bootstrap iterations used for S12: %d", bootstrap_iterations),
    "",
    "Generated manuscript-facing tables:",
    "- `TABLE_S6_agp_ibd_cleaner_control_sensitivity.tsv` / `.tex`",
    "- `TABLE_S7_ibd_dominance_lineage_retained_feature_profiles.tsv` / `.tex`",
    "- `TABLE_S8_ibd_antibiotic_status_sensitivity.tsv` / `.tex`",
    "- `TABLE_S9_absorb_policy_selected_ibd_robustness.tsv` / `.tex`",
    "- `TABLE_S10_selected_ibd_label_taxonomic_ranks.tsv` / `.tex`",
    "- `TABLE_S11A_rebuilt_external_validation_summary.tsv` / `.tex`",
    "- `TABLE_S11B_transfer_external_validation_summary.tsv` / `.tex`",
    "- `TABLE_S11C_external_validation_cohort_context.tsv` / `.tex`",
    "- `TABLE_S11D_external_ibd_cohort_characterization.tsv` / `.tex`",
    "- `TABLE_S12_external_subject_level_sensitivity.tsv` / `.tex`",
    "- `TABLE_S13_external_label_transfer_mapping_balance.tsv` / `.tex`",
    "- `TABLE_S13A_S13C_halfvarson_clinical_followup_summary.tsv` / `.tex`",
    "- `TABLE_S13D_halfvarson_clinical_followup_subject_sensitivity.tsv` / `.tex`",
    "- `TABLE_S13E_gevers_medication_followup.tsv` / `.tex`",
    "",
    sprintf("Cleaner-control rows: %d", nrow(cleaner_control)),
    sprintf("Retained-feature profile rows: %d", nrow(profiles)),
    sprintf("Antibiotic-status rows: %d", nrow(antibiotic)),
    sprintf("S11A rows: %d", nrow(overview$s11a)),
    sprintf("S11B rows: %d", nrow(overview$s11b)),
    sprintf("S11C rows: %d", nrow(overview$s11c)),
    sprintf("S11D rows: %d", nrow(overview$s11d)),
    sprintf("Subject-level bootstrap rows: %d", nrow(subject)),
    sprintf("Mapping-balance rows: %d", nrow(mapping)),
    sprintf("Halfvarson clinical follow-up rows: %d", nrow(clinical_followup$summary)),
    sprintf("Halfvarson clinical subject-sensitivity rows: %d", nrow(clinical_followup$subject_summary)),
    sprintf("Gevers medication follow-up rows: %d", nrow(gevers_medication)),
    sprintf("Absorb-policy rows: %d", nrow(policy)),
    sprintf("Taxonomic-rank rows: %d", nrow(ranks))
  )
  writeLines(lines, file.path(tables_dir, "TABLE_S8_S13_reviewer_sensitivity_summary.md"))
}

cleaner_control <- build_cleaner_control_sensitivity()
profiles <- build_retained_feature_profiles()
antibiotic <- build_antibiotic_sensitivity()
policy <- build_absorb_policy_sensitivity()
ranks <- build_taxonomic_rank_table()
overview <- build_external_validation_overview()
subject <- build_subject_level_diagnostics()
mapping <- build_mapping_balance()
clinical_followup <- build_halfvarson_clinical_followup()
gevers_medication <- build_gevers_medication_followup()
write_summary(cleaner_control, profiles, antibiotic, overview, subject, mapping, policy, ranks, clinical_followup, gevers_medication)

message("Wrote reviewer sensitivity assets to ", tables_dir)
