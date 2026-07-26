## Regenerate the AGP depth-1 DCST bar plot embedded in README.Rmd.

devtools::load_all(".", quiet = TRUE)
load("data/agp_gut.rda")

filtered <- filter.asv(
  agp_gut$counts,
  min.lib = 1000,
  prev.prop = 0.05,
  min.count = 2
)
normalized <- normalize.linf(filtered$counts)
dcsts <- linf.csts(normalized, n0 = 30)
sizes <- sort(table(dcsts$cell.label), decreasing = TRUE)

png(
  "man/figures/readme-dcst-barplot.png",
  width = 1400,
  height = 900,
  res = 150
)
par(mar = c(11, 5, 3, 1))
bars <- barplot(
  sizes,
  col = ifelse(names(sizes) == "RARE_DOMINANT", "#e74c3c", "#3498db"),
  las = 2,
  cex.names = 0.85,
  ylab = "Number of samples",
  main = "Gut DCST Size Distribution (depth 1)"
)
text(bars, sizes, labels = sizes, pos = 3, cex = 0.8)
dev.off()
