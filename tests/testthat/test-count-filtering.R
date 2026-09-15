test_that("filtering rejects invalid data and scalar thresholds", {
  X <- rbind(c(10, 1), c(10, 2))
  for (bad in c(-1, NA, NaN, Inf)) {
    Y <- X; Y[1, 2] <- bad
    expect_error(filter.asv(Y, verbose = FALSE), "finite and nonnegative")
  }
  expect_error(filter.asv(matrix("1", 2, 2)), "numeric")
  expect_error(filter.asv(data.frame(a = c("1", "2"))), "numeric")
  invalid <- list(min.lib = c(-1, 1.5, NA, Inf), prev.prop = c(0, 1.1, NA, Inf),
                  min.count = c(0, 1.5, NA, Inf), min.rel = c(0, 1, NA, Inf),
                  min.feat.total = c(-1, 1.5, NA, Inf))
  for (arg in names(invalid)) {
    for (value in invalid[[arg]]) {
      expect_error(do.call(filter.asv, c(list(S.counts = X, verbose = FALSE),
                                         setNames(list(value), arg))), arg, fixed = TRUE)
    }
    expect_error(do.call(filter.asv, c(list(S.counts = X, verbose = FALSE),
                                       setNames(list(c(1, 2)), arg))), arg, fixed = TRUE)
  }
})

test_that("relative prevalence handles zero-total rows explicitly", {
  X <- rbind(zero = c(A = 0, B = 0), positive = c(A = 10, B = 2))
  out <- filter.asv(X, min.lib = 0, min.rel = .1, verbose = FALSE)
  expect_equal(out$counts, X[2, , drop = FALSE])
  expect_equal(unname(out$kept.sample.idx), 2L)
  expect_equal(unname(out$prevalence), c(1, 1))
  expect_equal(unname(rowSums(out$rel)), 1)
  # Zero rows remain in the prevalence denominator until the final row filter.
  strict <- filter.asv(X, min.lib = 0, min.rel = .1, prev.prop = 1, verbose = FALSE)
  expect_equal(dim(strict$counts), c(0L, 0L))
  expect_equal(strict$thresholds$prev.thld, 2)
})

test_that("all filtering exits have coherent fields and indices", {
  X <- rbind(s1 = c(A = 5, B = 0), s2 = c(A = 0, B = 1), s3 = c(A = 4, B = 0))
  normal <- filter.asv(X, min.lib = 0, prev.prop = .5, min.count = 1, verbose = FALSE)
  expect_equal(normal$counts, X[c(1, 3), 1, drop = FALSE])
  expect_equal(normal$counts, X[normal$kept.sample.idx, normal$kept.feature.idx, drop = FALSE])
  expect_equal(unname(rowSums(normal$rel)), c(1, 1))
  expect_warning(empty <- filter.asv(X, min.lib = 100, verbose = FALSE), "No data left")
  expect_identical(names(empty), names(normal))
  expect_equal(dim(empty$counts), c(0L, 0L))
  expect_equal(empty$counts, X[empty$kept.sample.idx, empty$kept.feature.idx, drop = FALSE])
  expect_equal(empty$thresholds$prev.thld, 0)
  expect_length(empty$prevalence, 0)
  none <- filter.asv(X, min.lib = 0, min.count = 100, verbose = FALSE)
  expect_identical(names(none), names(normal))
  expect_equal(dim(none$rel), c(0L, 0L))
  fractional <- filter.asv(data.frame(A = c(1.5, 2.5)), min.lib = 0, min.count = 1, verbose = FALSE)
  expect_equal(unname(as.vector(fractional$counts)), c(1.5, 2.5))
  expect_equal(as.vector(fractional$rel), c(1, 1))
  one <- filter.asv(matrix(2, 1, 1), min.lib = 0, verbose = FALSE)
  expect_equal(dim(one$rel), c(1L, 1L))
})
