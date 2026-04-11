test_that("linf.landmarks computes feature-specific landmarks at depth 1", {
  M <- rbind(
    s1 = c(1.00, 0.20, 0.10),
    s2 = c(0.90, 0.40, 0.30),
    s3 = c(0.70, 0.60, 0.10),
    s4 = c(0.20, 1.00, 0.10),
    s5 = c(0.30, 0.90, 0.20)
  )
  ids <- c("asv_1", "asv_2", "asv_3")
  labels <- c("A", "B", "C")

  csts <- linf.csts(
    M,
    feature.ids = ids,
    feature.labels = labels,
    n0 = 2,
    low.freq.policy = "absorb"
  )
  out <- linf.landmarks(
    M,
    csts,
    landmark.types = c("endpoint.max", "endpoint.min", "mean.rep", "median.rep")
  )

  a.rows <- out$landmarks[out$landmarks$cell.id == "asv_1", , drop = FALSE]

  expect_s3_class(out, "linf.landmarks")
  expect_equal(a.rows$landmark.type, c("endpoint.max", "endpoint.min", "mean.rep", "median.rep"))
  expect_equal(a.rows$point.name, c("s1", "s3", "s2", "s2"))
  expect_true(all(a.rows$target.feature.id == "asv_1"))
  expect_true(all(a.rows$target.feature.label == "A"))
  expect_equal(a.rows$target.value[a.rows$landmark.type == "mean.rep"], mean(c(1.0, 0.9, 0.7)))
  expect_equal(a.rows$abs.deviation[a.rows$landmark.type == "mean.rep"], abs(0.9 - mean(c(1.0, 0.9, 0.7))))
})

test_that("linf.csts can attach landmarks and skip rare buckets", {
  A <- matrix(c(5, 1, 0), nrow = 6, ncol = 3, byrow = TRUE)
  B <- matrix(c(1, 5, 0), nrow = 2, ncol = 3, byrow = TRUE)
  C <- matrix(c(1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
  M <- rbind(A, B, C)
  M <- sweep(M, 1, rowSums(M), "/")
  colnames(M) <- c("Dom1", "Dom2", "Dom3")

  out <- linf.csts(
    M,
    n0 = 5,
    low.freq.policy = "pure",
    return.landmarks = TRUE,
    landmark.types = "endpoint.max",
    landmark.view = "rare"
  )

  rare.row <- out$landmarks$cells[out$landmarks$cells$cell.id == out$rare.label, , drop = FALSE]

  expect_s3_class(out$landmarks, "linf.landmarks")
  expect_identical(out$landmarks$view, "rare")
  expect_equal(nrow(rare.row), 1L)
  expect_false(rare.row$landmarks.computable)
  expect_true(all(out$landmarks$landmarks$cell.id != out$rare.label))
})

test_that("linf.landmarks uses leaf feature ids for refined depth-2 cells", {
  M <- rbind(
    r1 = c(1.00, 0.80, 0.20),
    r2 = c(0.95, 0.70, 0.30),
    r3 = c(0.90, 0.20, 0.70),
    r4 = c(0.85, 0.10, 0.80),
    r5 = c(0.20, 1.00, 0.10),
    r6 = c(0.30, 0.95, 0.10)
  )
  ids <- c("asv_1", "asv_2", "asv_3")
  labels <- c("A", "B", "C")

  d1 <- linf.csts(
    M,
    feature.ids = ids,
    feature.labels = labels,
    n0 = 2,
    low.freq.policy = "pure"
  )
  d2 <- refine.linf.csts(
    M,
    d1,
    n0 = 2,
    refinement.factor = 2,
    low.freq.policy = "pure",
    verbose = FALSE
  )
  out <- linf.landmarks(M, d2, depth = 2, view = "rare", landmark.types = "endpoint.max")

  c.row <- out$landmarks[out$landmarks$cell.id == "asv_1__asv_3", , drop = FALSE]

  expect_true("cst.id.levels.rare" %in% names(d2))
  expect_true("cst.id.levels.absorb" %in% names(d2))
  expect_equal(c.row$target.feature.id, "asv_3")
  expect_equal(c.row$target.feature.label, "C")
  expect_equal(c.row$point.name, "r4")
})

test_that("switching CST views updates ids and drops stale attached landmarks", {
  M <- rbind(
    s1 = c(1.0, 0.2, 0.1),
    s2 = c(0.9, 0.3, 0.1),
    s3 = c(0.1, 1.0, 0.1)
  )
  colnames(M) <- c("A", "B", "C")

  csts <- linf.csts(
    M,
    n0 = 3,
    low.freq.policy = "pure",
    return.landmarks = TRUE
  )
  collapsed <- collapse.rare(csts)

  expect_null(collapsed$landmarks)
  expect_identical(collapsed$cell.id, collapsed$cell.id.absorb)
  expect_identical(collapsed$cst.id.levels[[collapsed$cst.depth]], collapsed$cst.id.levels.absorb[[collapsed$cst.depth]])
})

test_that('low.freq.policy = "rare" is accepted as a deprecated alias for "pure"', {
  M <- rbind(
    s1 = c(1.0, 0.2, 0.1),
    s2 = c(0.9, 0.3, 0.1),
    s3 = c(0.1, 1.0, 0.1)
  )
  colnames(M) <- c("A", "B", "C")

  expect_warning(
    csts <- linf.csts(M, n0 = 2, low.freq.policy = "rare"),
    'deprecated; use "pure" instead'
  )

  expect_identical(csts$low.freq.policy, "pure")
  expect_identical(csts$cell.label, csts$cell.label.rare)
})
