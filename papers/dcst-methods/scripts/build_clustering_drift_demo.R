#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(cluster)
  library(ggplot2)
  library(scales)
})

root <- "/Users/pgajer/current_projects/linf"
paper_root <- file.path(root, "papers", "dcst-methods")
out_dir <- file.path(paper_root, "assets", "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(file.path(root, "data", "agp_gut.rda"))

counts <- agp_gut[["counts"]]
meta <- agp_gut[["meta"]]
rownames(meta) <- rownames(counts)

prop <- sweep(counts, 1, rowSums(counts), "/")
prop <- prop[, colSums(prop) > 0, drop = FALSE]

jsd_dist <- function(x) {
  n <- nrow(x)
  d <- matrix(0, n, n, dimnames = list(rownames(x), rownames(x)))
  entropy <- -rowSums(ifelse(x > 0, x * log(x), 0))
  for (i in seq_len(n - 1L)) {
    j <- seq.int(i + 1L, n)
    m <- sweep(x[j, , drop = FALSE], 2, x[i, ], "+") / 2
    mid_entropy <- -rowSums(ifelse(m > 0, m * log(m), 0))
    jsd <- pmax(mid_entropy - 0.5 * entropy[i] - 0.5 * entropy[j], 0)
    d[i, j] <- d[j, i] <- sqrt(jsd)
  }
  as.dist(d)
}

adjusted_rand <- function(a, b) {
  tab <- table(a, b)
  choose2 <- function(z) z * (z - 1) / 2
  n <- length(a)
  index <- sum(choose2(tab))
  row_index <- sum(choose2(rowSums(tab)))
  col_index <- sum(choose2(colSums(tab)))
  expected <- row_index * col_index / choose2(n)
  maximum <- 0.5 * (row_index + col_index)
  if (isTRUE(all.equal(maximum, expected))) {
    return(1)
  }
  (index - expected) / (maximum - expected)
}

permute <- function(x) {
  if (length(x) == 1L) {
    return(list(x))
  }
  unlist(
    lapply(seq_along(x), function(i) {
      lapply(permute(x[-i]), function(y) c(x[i], y))
    }),
    recursive = FALSE
  )
}

best_matched_agreement <- function(a, b) {
  if (identical(as.character(a), as.character(b))) {
    return(1)
  }
  tab <- table(factor(a), factor(b))
  k <- max(nrow(tab), ncol(tab))
  padded <- matrix(0, k, k)
  padded[seq_len(nrow(tab)), seq_len(ncol(tab))] <- tab
  best <- max(vapply(
    permute(seq_len(k)),
    function(p) sum(padded[cbind(seq_len(k), p)]),
    numeric(1)
  ))
  best / length(a)
}

set.seed(20260428)
target_ids <- sort(sample(rownames(prop), 120))
context_ids <- sort(sample(setdiff(rownames(prop), target_ids), 120))
combined_ids <- c(target_ids, context_ids)

target_prop <- prop[target_ids, , drop = FALSE]
combined_prop <- prop[combined_ids, , drop = FALSE]
target_dist <- jsd_dist(target_prop)
combined_dist <- jsd_dist(combined_prop)

pam_labels <- function(d, k) {
  as.integer(cluster::pam(d, k = k, diss = TRUE, variant = "faster")$clustering)
}

target_k3 <- pam_labels(target_dist, 3)
target_k4 <- pam_labels(target_dist, 4)
target_k5 <- pam_labels(target_dist, 5)
combined_k4 <- pam_labels(combined_dist, 4)
combined_k4_target <- combined_k4[match(target_ids, combined_ids)]

dcst_labels <- meta[target_ids, "dcst_depth2"]

comparisons <- data.frame(
  comparison = c(
    "PAM k=3 vs k=4\nsame 120 AGP samples",
    "PAM k=4 vs k=5\nsame 120 AGP samples",
    "PAM k=4\n120 samples alone vs +120 context",
    "Frozen dCST depth 2\nsame 120 AGP samples"
  ),
  ari = c(
    adjusted_rand(target_k3, target_k4),
    adjusted_rand(target_k4, target_k5),
    adjusted_rand(target_k4, combined_k4_target),
    adjusted_rand(dcst_labels, dcst_labels)
  ),
  matched_agreement = c(
    best_matched_agreement(target_k3, target_k4),
    best_matched_agreement(target_k4, target_k5),
    best_matched_agreement(target_k4, combined_k4_target),
    best_matched_agreement(dcst_labels, dcst_labels)
  ),
  stringsAsFactors = FALSE
)

comparisons$comparison <- factor(comparisons$comparison, levels = comparisons$comparison)

plot_df <- rbind(
  data.frame(comparison = comparisons$comparison, metric = "Adjusted Rand index", value = comparisons$ari),
  data.frame(comparison = comparisons$comparison, metric = "Best matched sample agreement", value = comparisons$matched_agreement)
)
plot_df$metric <- factor(plot_df$metric, levels = c("Adjusted Rand index", "Best matched sample agreement"))

write.table(
  comparisons,
  file.path(paper_root, "assets", "tables", "clustering_drift_demo_metrics.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE
)

p <- ggplot(plot_df, aes(x = value, y = comparison, fill = metric)) +
  geom_col(position = position_dodge(width = 0.72), width = 0.62, color = "grey20", linewidth = 0.25) +
  geom_text(
    aes(label = sprintf("%.2f", value)),
    position = position_dodge(width = 0.72),
    hjust = -0.10,
    size = 3.6,
    color = "grey15"
  ) +
  scale_x_continuous(limits = c(0, 1.08), breaks = seq(0, 1, by = 0.25), labels = label_number(accuracy = 0.01)) +
  scale_fill_manual(values = c("Adjusted Rand index" = "#4c78a8", "Best matched sample agreement" = "#f58518")) +
  labs(
    title = "Real-data illustration of clustering-label dependence",
    subtitle = "AGP example data: PAM labels vary with k/comparison set; frozen dCST labels are unchanged",
    x = "Label stability metric for target samples",
    y = NULL,
    fill = NULL
  ) +
  coord_cartesian(clip = "off") +
  theme_minimal(base_size = 11) +
  theme(
    plot.title = element_text(face = "bold", size = 14, margin = margin(b = 4)),
    plot.subtitle = element_text(size = 10.5, color = "grey25", margin = margin(b = 12)),
    axis.text.y = element_text(size = 9.4, color = "grey15"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "top",
    legend.justification = "left",
    plot.margin = margin(12, 28, 12, 12)
  )

png_out <- file.path(out_dir, "FIGURE_S2_clustering_drift_demo.png")
pdf_out <- file.path(out_dir, "FIGURE_S2_clustering_drift_demo.pdf")
ggsave(png_out, p, width = 8.8, height = 5.3, dpi = 320)
pdf(pdf_out, width = 8.8, height = 5.3, version = "1.4", useDingbats = FALSE)
print(p)
dev.off()

message("Wrote ", png_out)
message("Wrote ", pdf_out)
