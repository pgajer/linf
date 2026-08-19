test_that("public dCST objects use dominance-lineage terminology", {
  M <- rbind(
    s1 = c(A = 1.0, B = 0.2),
    s2 = c(A = 0.9, B = 0.3),
    s3 = c(A = 0.1, B = 1.0)
  )

  fit <- linf.csts(M, n0 = 2, low.freq.policy = "pure")

  expect_true(all(c(
    "depth1.feature.index", "depth1.feature.id", "depth1.feature.label",
    "lineage.id", "lineage.label", "lineage.ids", "lineage.labels", "depth"
  ) %in% names(fit)))
  expect_false(any(c(
    "cell.index", "cell.id", "cell.label", "kept.cells.idx",
    "cst.id.levels", "cst.levels", "cst.depth"
  ) %in% names(fit)))

  landmarks <- linf.landmarks(M, fit, view = "pure")
  expect_true("lineages" %in% names(landmarks))
  expect_false("cells" %in% names(landmarks))
  expect_true(all(c("lineage.id", "lineage.label", "lineage.size") %in%
                    names(landmarks$lineages)))
})

test_that("removed public aliases are not present", {
  ns <- asNamespace("linf")

  expect_false(exists("linf.cells", envir = ns, inherits = FALSE))
  expect_false(exists("collapse.rare", envir = ns, inherits = FALSE))
  expect_false(exists("expand.rare", envir = ns, inherits = FALSE))
})
