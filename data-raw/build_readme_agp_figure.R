## Regenerate the AGP absorb-policy depth-1 dCST bar plot embedded in README.Rmd.

devtools::load_all(".", quiet = TRUE)
load("data/agp_gut.rda")

filtered <- filter.asv(
  agp_gut$counts,
  min.lib = 1000,
  prev.prop = 0.05,
  min.count = 2
)
normalized <- normalize.linf(filtered$counts)
dcsts <- linf.csts(
  normalized,
  n0 = 30,
  low.freq.policy = "absorb",
  return.diagnostics = TRUE
)
sizes <- sort(table(dcsts$lineage.label), decreasing = TRUE)
stopifnot(
  sum(sizes) == nrow(normalized),
  !"RARE_DOMINANT" %in% names(sizes)
)

short_taxon_label <- function(x) {
  ranks <- strsplit(x, ";", fixed = TRUE)[[1L]]
  labels <- sub("^[a-z]__", "", ranks)
  labels <- labels[nzchar(labels) & labels != "__"]
  if (!length(labels)) x else labels[[length(labels)]]
}

names(sizes) <- make.unique(
  vapply(names(sizes), short_taxon_label, character(1L))
)

png(
  "man/figures/readme-dcst-barplot.png",
  width = 1400,
  height = 900,
  res = 150
)
par(mar = c(8, 5, 4, 1))
bars <- barplot(
  sizes,
  col = "#3498db",
  border = "#1f618d",
  las = 2,
  cex.names = 0.9,
  ylim = c(0, max(sizes) * 1.15),
  ylab = "Number of samples",
  main = "Gut dCST Size Distribution (absorb policy, depth 1)"
)
text(bars, sizes, labels = sizes, pos = 3, cex = 0.8)
dev.off()
