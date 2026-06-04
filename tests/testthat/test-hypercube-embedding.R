test_that("linf.hypercube.embedding handles finite ratios and boundary rows", {
  X <- rbind(
    finite = c(A = 2, B = 1, C = 1),
    boundary = c(A = 0, B = 2, C = 1),
    zero = c(A = 0, B = 0, C = 0),
    reference_only = c(A = 1, B = 0, C = 0)
  )

  emb <- linf.hypercube.embedding(X, reference = "A", lambda = log(2))

  expect_equal(dim(emb), c(4L, 2L))
  expect_equal(colnames(emb), c("B_rel_A", "C_rel_A"))
  expect_equal(unname(emb["finite", ]), c(0.5, 0.5), tolerance = 1e-12)
  expect_equal(unname(emb["boundary", ]), c(1, 0.5), tolerance = 1e-12)
  expect_equal(unname(emb["zero", ]), c(0, 0), tolerance = 1e-12)
  expect_equal(unname(emb["reference_only", ]), c(0, 0), tolerance = 1e-12)

  expect_equal(attr(emb, "reference.index"), 1L)
  expect_equal(attr(emb, "reference.id"), "A")
  expect_equal(attr(emb, "lambda"), log(2))
  expect_equal(attr(emb, "lambda.policy"), "fixed")
  expect_equal(attr(emb, "finite.reference.count"), 2L)
  expect_equal(attr(emb, "zero.reference.count"), 2L)
})

test_that("linf.hypercube.embedding is projectively scale invariant", {
  X <- rbind(
    c(A = 3, B = 6, C = 9),
    c(A = 0, B = 4, C = 8)
  )
  Y <- X * c(2, 7)

  emb.x <- linf.hypercube.embedding(X, reference = "A", lambda = 0.25)
  emb.y <- linf.hypercube.embedding(Y, reference = "A", lambda = 0.25)

  expect_equal(unname(emb.x), unname(emb.y), tolerance = 1e-12)
})

test_that("linf.hypercube.embedding resolves feature labels and auto lambda", {
  X <- rbind(
    c(2, 1, 1),
    c(4, 2, 2)
  )
  colnames(X) <- c("asv1", "asv2", "asv3")

  emb <- linf.hypercube.embedding(
    X,
    reference = "Li",
    feature.labels = c("Li", "Lc", "Gv"),
    sigma.quantile = 1,
    sigma.target = 0.5
  )

  expect_equal(attr(emb, "reference.id"), "asv1")
  expect_equal(attr(emb, "reference.label"), "Li")
  expect_equal(attr(emb, "lambda"), log(2), tolerance = 1e-12)
  expect_equal(unname(emb[1, ]), c(0.5, 0.5), tolerance = 1e-12)
})

test_that("linf.hypercube.embedding agrees for dense and sparse inputs", {
  skip_if_not_installed("Matrix")

  X <- rbind(
    c(A = 2, B = 1, C = 1),
    c(A = 0, B = 2, C = 1),
    c(A = 1, B = 0, C = 0)
  )

  dense <- linf.hypercube.embedding(X, reference = "A", lambda = log(2))
  sparse <- linf.hypercube.embedding(Matrix::Matrix(X, sparse = TRUE),
                                     reference = "A", lambda = log(2))

  expect_equal(unname(dense), unname(sparse), tolerance = 1e-12)
  expect_equal(colnames(sparse), colnames(dense))
})

test_that("linf.hypercube.embedding validates inputs", {
  expect_error(
    linf.hypercube.embedding(matrix(c(1, -1), nrow = 1), reference = 1),
    "negative entries"
  )
  expect_error(
    linf.hypercube.embedding(matrix(c(1, 2), nrow = 1), reference = 3),
    "out of range"
  )
  expect_error(
    linf.hypercube.embedding(matrix(c(1, 2), nrow = 1), reference = 1, lambda = 0),
    "lambda"
  )
})
