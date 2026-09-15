test_that("matrix consumers check both dimensions and named fitted axes", {
  X <- rbind(s1 = c(A = 10, B = 1), s2 = c(A = 9, B = 2),
             s3 = c(A = 1, B = 10), s4 = c(A = 2, B = 9))
  for (sparse in c(FALSE, TRUE)) {
    M <- if (sparse) Matrix::Matrix(X, sparse = TRUE) else X
    f <- linf.csts(M, feature.ids = c("id_A", "id_B"), n0 = 2)
    expect_identical(f$input.dimnames, dimnames(X))
    for (consume in list(
      function(M, f) linf.landmarks(M, f),
      function(M, f) refine.linf.csts(M, f, lineages.to.refine = "id_A", n0 = 2, verbose = FALSE)
    )) {
      expect_no_error(consume(M, f))
      expect_error(consume(M[, 2:1], f), "column names/order.*align")
      expect_error(consume(M[4:1, ], f), "row names/order.*align")
      expect_error(consume(M[1:2, ], f), "row count.*expected 4")
      expect_error(consume(M[, 1, drop = FALSE], f), "column count.*expected 2")
      expect_no_error(consume(unname(M), f))
      legacy <- f; legacy$input.dimnames <- NULL
      expect_no_error(consume(M, legacy))
      expect_error(consume(M[1:2, ], legacy), "row count")
    }
  }
})

test_that("only unique nonempty query matching keys can drive transfer", {
  X <- rbind(c(A = 10, B = 1), c(A = 1, B = 10))
  f <- linf.csts(X, n0 = 1)
  Q <- matrix(c(0, 9, 2), nrow = 1)
  for (keys in list(c("A", "A", "B"), c("A", NA, "B"), c("A", "", "B"))) {
    expect_error(transfer.dcsts(Q, f, feature.ids = keys), "keys|unique")
    expect_error(transfer.dcsts(Q, f, feature.labels = keys, match.by = "feature.labels"), "keys|unique")
  }
  expect_no_error(transfer.dcsts(Q, f, feature.ids = c("A", "C", "B"),
                                feature.labels = c("display", "display", "display")))
  expect_equal(transfer.dcsts(X[, 2:1], f)$assignment.ids,
               transfer.dcsts(X, f)$assignment.ids)
})
