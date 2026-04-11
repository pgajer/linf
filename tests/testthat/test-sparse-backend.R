test_that("normalize.linf returns the same values under dense and sparse backends", {
  M <- rbind(
    s1 = c(5, 2, 0, 0),
    s2 = c(0, 0, 0, 0),
    s3 = c(1, 4, 4, 0),
    s4 = c(0, 3, 0, 6)
  )

  dense <- normalize.linf(M, backend = "dense")
  sparse <- normalize.linf(Matrix::Matrix(M, sparse = TRUE), backend = "auto")

  expect_true(inherits(sparse, "dgCMatrix"))
  expect_equal(as.matrix(sparse), dense)
})

test_that("linf.cells matches between dense and sparse backends", {
  M <- rbind(
    s1 = c(1, 0.5, 0),
    s2 = c(0, 0, 0),
    s3 = c(0.2, 0.8, 0.8),
    s4 = c(0, 0.4, 0.9)
  )
  ids <- c("asv_1", "asv_2", "asv_3")
  labels <- c("A", "B", "C")

  dense <- linf.cells(M, feature.ids = ids, feature.labels = labels, backend = "dense")
  sparse <- linf.cells(Matrix::Matrix(M, sparse = TRUE), feature.ids = ids, feature.labels = labels, backend = "auto")

  expect_equal(sparse$index, dense$index)
  expect_equal(sparse$id, dense$id)
  expect_equal(sparse$label, dense$label)
  expect_equal(sparse$observed.id.levels, dense$observed.id.levels)
})

test_that("linf.csts and refinement agree across dense and sparse backends", {
  M <- rbind(
    r1 = c(1.00, 0.70, 0.10, 0.00),
    r2 = c(0.95, 0.75, 0.10, 0.00),
    r3 = c(0.90, 0.20, 0.80, 0.00),
    r4 = c(0.88, 0.10, 0.82, 0.00),
    r5 = c(0.10, 1.00, 0.05, 0.00),
    r6 = c(0.12, 0.95, 0.08, 0.00),
    r7 = c(0.05, 0.10, 0.05, 1.00),
    r8 = c(0.06, 0.12, 0.04, 0.95)
  )
  ids <- paste0("asv_", seq_len(ncol(M)))
  labels <- LETTERS[seq_len(ncol(M))]

  dense.d1 <- linf.csts(
    M,
    feature.ids = ids,
    feature.labels = labels,
    n0 = 2,
    low.freq.policy = "absorb",
    backend = "dense"
  )
  sparse.input <- Matrix::Matrix(M, sparse = TRUE)
  sparse.d1 <- linf.csts(
    sparse.input,
    feature.ids = ids,
    feature.labels = labels,
    n0 = 2,
    low.freq.policy = "absorb",
    backend = "auto"
  )

  expect_identical(sparse.d1$matrix.backend, "sparse")
  expect_equal(sparse.d1$cell.id, dense.d1$cell.id)
  expect_equal(sparse.d1$cell.label, dense.d1$cell.label)
  expect_equal(sparse.d1$cell.label.absorb, dense.d1$cell.label.absorb)
  expect_equal(sparse.d1$cell.label.rare, dense.d1$cell.label.rare)

  dense.d2 <- refine.linf.csts(
    M,
    dense.d1,
    n0 = 2,
    refinement.factor = 2,
    low.freq.policy = "absorb",
    verbose = FALSE,
    backend = "dense"
  )
  sparse.d2 <- refine.linf.csts(
    sparse.input,
    sparse.d1,
    n0 = 2,
    refinement.factor = 2,
    low.freq.policy = "absorb",
    verbose = FALSE
  )

  expect_identical(sparse.d2$matrix.backend, "sparse")
  expect_equal(sparse.d2$cst.levels[[2]], dense.d2$cst.levels[[2]])
  expect_equal(sparse.d2$cst.id.levels[[2]], dense.d2$cst.id.levels[[2]])
})

test_that("linf.landmarks and the landmark pipeline accept sparse matrices", {
  M <- rbind(
    s1 = c(10, 8, 1),
    s2 = c(9, 7, 2),
    s3 = c(8, 2, 7),
    s4 = c(7, 1, 8),
    s5 = c(1, 10, 2),
    s6 = c(2, 9, 1)
  )
  ids <- c("asv_1", "asv_2", "asv_3")
  labels <- c("A", "B", "C")
  sparse.M <- Matrix::Matrix(M, sparse = TRUE)

  sparse.rel <- normalize.linf(sparse.M, backend = "auto")
  sparse.d1 <- linf.csts(
    sparse.rel,
    feature.ids = ids,
    feature.labels = labels,
    n0 = 2,
    low.freq.policy = "pure",
    return.landmarks = TRUE,
    backend = "auto"
  )

  expect_s3_class(sparse.d1$landmarks, "linf.landmarks")

  pipe.out <- linf.dcst.landmark.pipeline(
    sparse.M,
    feature.ids = ids,
    feature.labels = labels,
    n0.depth1 = 2,
    n0.depth2 = 2,
    refinement.factor = 2,
    low.freq.policy = "pure",
    landmark.view = "absorb",
    verbose = FALSE,
    backend = "auto"
  )

  expect_identical(pipe.out$dcst.depth1$matrix.backend, "sparse")
  expect_s3_class(pipe.out$landmarks.depth2, "linf.landmarks")
})
