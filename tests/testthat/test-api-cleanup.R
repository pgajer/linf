test_that("removed 0.3.0 functions and helpers are absent", {
  ns <- asNamespace("linf")

  removed <- c(
    "asv.to.linf.csts",
    "latex.linf.csts",
    "refine.linf.csts.iter",
    "escape.latex",
    "df.to.latex"
  )
  expect_false(any(vapply(
    removed,
    exists,
    logical(1),
    envir = ns,
    inherits = FALSE
  )))
})

test_that("refine.linf.csts supports explicit and repeated refinement", {
  M <- rbind(
    s1 = c(A = 1.0, B = 0.8, C = 0.1),
    s2 = c(A = 1.0, B = 0.7, C = 0.2),
    s3 = c(A = 1.0, B = 0.6, C = 0.3),
    s4 = c(A = 0.2, B = 1.0, C = 0.8),
    s5 = c(A = 0.1, B = 1.0, C = 0.7),
    s6 = c(A = 0.3, B = 1.0, C = 0.6)
  )
  d1 <- linf.csts(M, n0 = 2, low.freq.policy = "absorb")

  d2 <- refine.linf.csts(
    M,
    d1,
    lineages.to.refine = "A",
    n0 = 2,
    low.freq.policy = "absorb",
    verbose = FALSE
  )
  expect_identical(d2$depth, 2L)
  expect_true(all(d2$lineage.id[1:3] == "A__B"))
  expect_true(all(d2$lineage.id[4:6] == "B"))

  d3 <- refine.linf.csts(
    M,
    d2,
    lineages.to.refine = "A__B",
    n0 = 2,
    low.freq.policy = "absorb",
    verbose = FALSE
  )
  expect_identical(d3$depth, 3L)
  expect_true(all(d3$lineage.id[1:3] == "A__B__C"))
  expect_true(all(d3$lineage.id[4:6] == "B"))

  expect_error(
    refine.linf.csts(
      M,
      d1,
      lineages.to.refine = "not-a-lineage",
      n0 = 2,
      verbose = FALSE
    ),
    "unknown lineage ID"
  )
})

test_that("pre-0.2 flat dCST objects are rejected", {
  M <- rbind(
    s1 = c(A = 1.0, B = 0.2),
    s2 = c(A = 0.9, B = 0.3),
    s3 = c(A = 0.1, B = 1.0)
  )
  fit <- linf.csts(M, n0 = 1)
  fit$lineage.labels <- NULL

  expect_error(
    print(fit),
    "incompatible dCST object"
  )
})
