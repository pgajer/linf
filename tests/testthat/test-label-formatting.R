test_that("linf.feature.labels applies abbreviations and global indices", {
  ids <- c("asv_1", "asv_4", "asv_5", "asv_6")
  tax <- c(
    "Lactobacillus iners",
    "Lactobacillus iners",
    "Megasphaera lornae",
    "Ca_Lachnocurva_vaginae"
  )
  abbr <- c(
    Lactobacillus = "L.",
    Streptococcus = "Strep",
    Gardnerella = "Gard.",
    Fannyhessea = "F.",
    Megasphaera = "Mega.",
    Ureaplasma = "U."
  )
  aliases <- c(
    "Ca. Lachnocurva vaginae" = "BVAB1",
    "Ca_Lachnocurva_vaginae" = "BVAB1"
  )

  out <- linf.feature.labels(
    feature.ids = ids,
    taxonomy = tax,
    abbreviations = abbr,
    aliases = aliases,
    duplicate.index = "global"
  )

  expect_equal(out, c("L. iners 1", "L. iners 4", "Mega. lornae", "BVAB1"))
})

test_that("linf.cells keeps ids separate from labels", {
  S <- rbind(c(0.7, 0.3), c(0, 0), c(0.1, 0.9))
  ids <- c("asv_1", "asv_4")
  labels <- c("L. iners 1", "L. iners 4")

  out <- linf.cells(S, feature.ids = ids, feature.labels = labels)

  expect_equal(out$index, c(1L, NA_integer_, 2L))
  expect_equal(out$id, c("asv_1", NA_character_, "asv_4"))
  expect_equal(out$label, c("L. iners 1", NA_character_, "L. iners 4"))
  expect_equal(out$id.levels, ids)
  expect_equal(out$levels, labels)
})

test_that("refine.linf.csts tracks id paths separately from label paths", {
  M <- rbind(
    c(0.90, 0.08, 0.02),
    c(0.88, 0.09, 0.03),
    c(0.87, 0.10, 0.03),
    c(0.86, 0.11, 0.03),
    c(0.05, 0.92, 0.03),
    c(0.06, 0.91, 0.03)
  )
  ids <- c("asv_1", "asv_4", "asv_5")
  labels <- c("L. iners 1", "L. iners 4", "Mega. lornae")

  d1 <- linf.csts(
    M,
    feature.ids = ids,
    feature.labels = labels,
    n0 = 2,
    low.freq.policy = "rare"
  )
  d2 <- refine.linf.csts(
    M,
    d1,
    n0 = 1,
    refinement.factor = 2,
    sep = "__",
    low.freq.policy = "rare",
    verbose = FALSE
  )

  expect_true("cst.id.levels" %in% names(d2))
  expect_true(all(grepl("^asv_1(__|$)|^asv_4(__|$)", d2$cst.id.levels[[2]][!is.na(d2$cst.id.levels[[2]])])))
  expect_true(all(grepl("^L\\. iners 1(__|$)|^L\\. iners 4(__|$)", d2$cst.levels[[2]][!is.na(d2$cst.levels[[2]])])))
})
