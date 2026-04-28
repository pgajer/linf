test_that("transfer.dcsts reproduces fitted labels on aligned training data", {
  X <- rbind(
    s1 = c(A = 10, B = 5, C = 1, D = 0),
    s2 = c(A = 9, B = 4, C = 2, D = 0),
    s3 = c(A = 8, B = 1, C = 6, D = 0),
    s4 = c(A = 7, B = 1, C = 5, D = 0),
    s5 = c(C = 10, A = 5, B = 1, D = 0),
    s6 = c(C = 9, A = 4, B = 1, D = 0),
    s7 = c(C = 8, A = 1, B = 6, D = 0),
    s8 = c(C = 7, A = 1, B = 5, D = 0)
  )
  X <- X[, c("A", "B", "C", "D")]
  M <- normalize.linf(X)
  fit <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")
  fit <- refine.linf.csts(M, fit, n0 = 2, low.freq.policy = "absorb", verbose = FALSE)

  tr <- transfer.dcsts(X, fit, match.by = "feature.labels")

  expect_equal(tr$assignment[, "depth1"], fit$cst.levels.absorb[[1]])
  expect_equal(tr$assignment[, "depth2"], fit$cst.levels.absorb[[2]])
})

test_that("transfer.dcsts aligns by feature metadata rather than column order", {
  X <- rbind(
    s1 = c(A = 10, B = 1, C = 0),
    s2 = c(A = 9, B = 2, C = 0),
    s3 = c(A = 1, B = 10, C = 0),
    s4 = c(A = 2, B = 9, C = 0)
  )
  fit <- linf.csts(normalize.linf(X), n0 = 2, low.freq.policy = "absorb")
  shuffled <- X[, c("C", "B", "A")]

  tr <- transfer.dcsts(shuffled, fit, match.by = "feature.labels")

  expect_equal(tr$assignment[, "depth1"], fit$cst.levels.absorb[[1]])
})

test_that("transfer.dcsts can use support tie-breaking in frozen transfer", {
  X <- rbind(
    a1 = c(A = 10, B = 1),
    a2 = c(A = 9, B = 1),
    a3 = c(A = 8, B = 1),
    b1 = c(A = 1, B = 10)
  )
  fit <- linf.csts(normalize.linf(X), n0 = 1, low.freq.policy = "absorb")
  tied <- rbind(new = c(A = 5, B = 5))

  tr <- transfer.dcsts(tied, fit, match.by = "feature.labels", tie.method = "support")

  expect_equal(unname(tr$assignment[1, "depth1"]), "A")
})
