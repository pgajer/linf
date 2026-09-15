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

test_that("finite extreme ratios keep their analytic direction and radius", {
  X <- rbind(overflow = c(A = 1e-320, B = 1, C = 0.5),
             tiny = c(A = 1, B = 1e-20, C = 0.5e-20),
             product = c(A = 1e200, B = 1e-200, C = 0.5e-200))
  out <- linf.hypercube.embedding(X, "A", lambda = 1)
  expect_equal(unname(out[1, ]), c(1, 0.5))
  expect_equal(unname(out[2, ]) / 1e-20, c(1.5, 0.75), tolerance = 1e-12)
  large.lambda <- linf.hypercube.embedding(X[3, , drop = FALSE], "A", lambda = 1e300)
  expect_equal(as.numeric(large.lambda) / 1e-100, c(1.5, 0.75), tolerance = 1e-12)
  expect_true(all(is.finite(out) & out >= 0 & out <= 1))
  sparse <- linf.hypercube.embedding(Matrix::Matrix(X, sparse = TRUE), "A", lambda = 1)
  expect_equal(sparse, out)
  expect_error(linf.hypercube.embedding(matrix(1, 2, 1), 1), "at least two features")
  expect_error(linf.hypercube.embedding(X, Inf), "out of range")
})

test_that("automatic scaling includes overflowed norms and is reusable in logs", {
  X <- rbind(c(A = 1e-320, B = 1e308, C = 0.5e308),
             c(A = 1e-320, B = 1e308, C = 0.5e308))
  out <- linf.hypercube.embedding(X, "A", sigma.target = 0.5)
  expect_equal(unname(out[1, ]), c(0.5, 0.25), tolerance = 1e-12)
  expect_true(is.na(attr(out, "lambda")))
  expect_true(is.finite(attr(out, "log.lambda")))
  fixed <- linf.hypercube.embedding(X, "A", log.lambda = attr(out, "log.lambda"))
  expect_equal(as.numeric(fixed), as.numeric(out))
  small <- linf.hypercube.embedding(X[, 3:1], "B", sigma.target = 0.5)
  expect_true(all(is.finite(small)))
  Y <- rbind(c(A = 1, B = 1, C = 1), c(A = 1, B = 4, C = 4))
  ordinary <- linf.hypercube.embedding(Y, "A", sigma.quantile = 0.25, sigma.target = 0.5)
  expect_equal(attr(ordinary, "lambda"), log(2) / 3.5, tolerance = 1e-12)
  expect_equal(linf.hypercube.embedding(Y * c(1e200, 1e-200), "A"),
               linf.hypercube.embedding(Y, "A"), tolerance = 1e-12)
  expect_error(linf.hypercube.embedding(Y, "A", lambda = 1, log.lambda = 0), "only one")
  expect_error(linf.hypercube.embedding(Y, "A", log.lambda = Inf), "finite")
})

test_that("tolerance and the zero-reference extension retain their meanings", {
  X <- rbind(c(A = 2, B = 0.2, C = 0), c(A = 0.1, B = 1, C = 0.5),
             c(A = 1, B = 0, C = 0), c(A = 0, B = 0, C = 0))
  out <- linf.hypercube.embedding(X, "A", lambda = 1, tol = 0.1)
  expect_equal(unname(out[1, ]), c(0, 0))
  expect_equal(unname(out[2, ]), c(1, 0.5))
  expect_true(all(out[3:4, ] == 0))
  approaching <- linf.hypercube.embedding(cbind(A = c(1, 0.1, 1e-320), B = 1, C = 0.5),
                                           "A", lambda = 1)
  expect_true(all(diff(approaching[, 1]) >= 0))
  expect_equal(unname(approaching[3, ]), c(1, 0.5))
})
