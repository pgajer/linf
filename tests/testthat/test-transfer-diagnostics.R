test_that("transfer IDs are independent of display labels and depth selection", {
  X <- rbind(s1 = c(A = 10, B = 5, C = 1), s2 = c(A = 9, B = 4, C = 1),
             s3 = c(A = 1, B = 0, C = 10), s4 = c(A = 1, B = 0, C = 9))
  fit <- function(labels) refine.linf.csts(X, linf.csts(X, n0 = 2, feature.labels = labels),
                                          lineages.to.refine = "A", n0 = 2, verbose = FALSE)
  f <- fit(c("Taxon A", "Taxon B", "Taxon C"))
  a <- transfer.dcsts(X, f, depth = c(2, 1))
  b <- transfer.dcsts(X, fit(c("Alpha", "Beta", "Gamma")), depth = c(2, 1))
  expect_identical(a$assignment.ids, b$assignment.ids)
  expect_false(identical(a$assignment, b$assignment))
  expect_equal(a$assignment.ids, a$all.depths.ids[, c(2, 1)])
  expect_equal(a$diagnostics$reason, rep("complete", 4))
  expect_equal(a$diagnostics$assigned.depth, rep(2L, 4))
  expect_true(all(is.na(a$diagnostics$stop.depth)))
  expect_identical(a$feature.match$matched, rep(TRUE, 3))
})

test_that("transfer diagnostics distinguish missing features from zero evidence", {
  X <- rbind(c(A = 10, B = 1), c(A = 9, B = 2))
  f <- refine.linf.csts(X, linf.csts(X, n0 = 2),
                        lineages.to.refine = "A", n0 = 2, verbose = FALSE)
  for (sparse in c(FALSE, TRUE)) {
    query <- function(x) if (sparse) Matrix::Matrix(x, sparse = TRUE) else x
    no.match <- transfer.dcsts(query(matrix(1, 1, 1, dimnames = list("q", "other"))), f)
    zero <- transfer.dcsts(query(matrix(c(0, 0), 1, dimnames = list("q", c("A", "B")))), f)
    missing.child <- transfer.dcsts(query(matrix(1, 1, 1, dimnames = list("q", "A"))), f)
    zero.child <- transfer.dcsts(query(matrix(c(1, 0), 1, dimnames = list("q", c("A", "B")))), f)
    expect_equal(no.match$diagnostics$reason, "missing_candidate_features")
    expect_equal(no.match$diagnostics$missing.candidates, 1L)
    expect_equal(zero$diagnostics$reason, "no_positive_candidate_values")
    expect_equal(missing.child$diagnostics$stop.depth, 2L)
    expect_equal(missing.child$diagnostics$reason, "missing_candidate_features")
    expect_equal(zero.child$diagnostics$reason, "no_positive_candidate_values")
    expect_equal(zero.child$diagnostics$assigned.depth, 1L)
    expect_equal(missing.child$feature.match$query.index, c(1L, NA_integer_))
  }
  tied <- transfer.dcsts(matrix(c(1, 1), 1), linf.csts(diag(2), n0 = 1))
  expect_equal(tied$diagnostics$n.tied.depths, 1L)
  empty <- linf.csts(diag(2), n0 = 3)
  expect_equal(transfer.dcsts(diag(2), empty)$diagnostics$reason, rep("no_candidates", 2))
  expect_equal(transfer.dcsts(diag(2), empty, view = "pure")$diagnostics$reason,
               rep("synthetic_only", 2))
})

test_that("history records each fit and refinement without inventing legacy settings", {
  X <- rbind(c(A = 10, B = 1), c(A = 9, B = 2))
  f <- linf.csts(X, n0 = 2, low.freq.policy = "absorb", tie.method = "random")
  g <- refine.linf.csts(X, f, n0 = 1, lineages.to.refine = "A", verbose = FALSE)
  expect_equal(g$n0, 2L)
  expect_equal(g$history[[1]]$tie.method, "random")
  expect_equal(g$history[[2]]$tie.method, "first")
  expect_equal(g$history[[2]]$n0, 1L)
  expect_equal(g$history[[2]]$source.view, "absorb")
  expect_equal(g$history[[2]]$low.freq.policy, "pure")
  expect_equal(g$history[[2]]$selected.lineage.ids, "A")
  expect_equal(g$history[[2]]$refined.node.ids, "1")
  expect_identical(dcst.view(g, "absorb")$history, g$history)
  tmp <- tempfile(); on.exit(unlink(tmp)); saveRDS(g, tmp)
  expect_identical(readRDS(tmp)$history, g$history)
  f$history <- NULL
  legacy <- refine.linf.csts(X, f, n0 = 3, verbose = FALSE)
  expect_identical(legacy$history[[1]], list(depth = 1L, available = FALSE))
  expect_identical(legacy$history[[2]]$refined.node.ids, character())
  expect_identical(legacy$history[[2]]$selection.mode, "automatic")
})
