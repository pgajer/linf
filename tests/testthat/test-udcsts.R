udcst.example <- function() {
  rbind(c(A = 9, B = 5, C = 1), c(A = 8, B = 4, C = 1),
        c(A = 4, B = 9, C = 1), c(A = 3, B = 8, C = 1),
        c(A = 1, B = 3, C = 9))
}

test_that("unordered supports pool opposite orders before filtering", {
  X <- udcst.example()
  fit <- linf.udcsts(X, level = 2, n0 = 3, low.freq.policy = "filter")
  expect_equal(fit$pure.id, c(rep("1,2", 4), "2,3"))
  expect_equal(fit$assignment, c(rep("1,2", 4), NA_character_))
  expect_identical(fit$keep, c(rep(TRUE, 4), FALSE))
  expect_equal(fit$status, c(rep("retained", 4), "filtered"))
  expect_equal(fit$feature.sets, rbind(c(1L, 2L), c(1L, 2L), c(1L, 2L),
                                     c(1L, 2L), c(2L, 3L)))
  expect_equal(fit$ordered.features[3, ], c(2L, 1L))
  expect_identical(X, udcst.example())
})

test_that("rare, filtered, and absorbed samples remain distinguishable", {
  X <- udcst.example()
  pure <- linf.udcsts(X, n0 = 3, low.freq.policy = "pure")
  absorb <- linf.udcsts(X, n0 = 3, low.freq.policy = "absorb")
  expect_true(all(pure$keep))
  expect_equal(pure$status[5], "rare")
  expect_true(is.na(pure$assignment[5]))
  expect_equal(absorb$status[5], "absorbed")
  expect_equal(absorb$assignment, rep("1,2", 5))
  expect_identical(absorb$support, pure$support)
  for (policy in c("pure", "filter")) {
    empty <- linf.udcsts(X, n0 = 6, low.freq.policy = policy)
    expect_length(empty$retained.ids, 0)
    expect_true(all(is.na(empty$assignment)))
    expect_identical(empty$keep, rep(policy == "pure", nrow(X)))
  }
  expect_error(linf.udcsts(X, n0 = 6, low.freq.policy = "absorb"),
               "no retained state")
})

test_that("support tiers have identical cumulative percentages", {
  X <- rbind(c(3, 0, 0), c(4, 0, 0), c(0, 3, 0), c(0, 4, 0), c(0, 0, 3))
  fit <- linf.udcsts(X, level = 1, n0 = 2, low.freq.policy = "filter")
  expect_identical(fit$frequency.table,
                   data.frame(ID = c("1", "2", "3"), freq = c(2L, 2L, 1L),
                              cumm.pct = c(80, 80, 100)))
  dcst <- linf.csts(X, n0 = 1, return.landmarks = FALSE)
  expect_equal(as.integer(fit$pure.id), dcst$provisional.feature.index)
})

test_that("ties and zero-valued ranks follow original columns", {
  X <- rbind(c(2, 2, 2, 1), c(0, 0, 4, 0))
  fit <- linf.udcsts(X, level = 3, n0 = 1)
  expect_equal(fit$ordered.features, rbind(c(1L, 2L, 3L), c(3L, 1L, 2L)))
  expect_equal(fit$pure.id, rep("1,2,3", 2))
  expect_identical(fit$positive.features, c(4L, 1L))
  expect_error(linf.udcsts(rbind(X, 0), n0 = 1), "all-zero rows")
  set.seed(712)
  before <- .Random.seed
  invisible(linf.udcsts(X, n0 = 1))
  expect_identical(.Random.seed, before)
})

test_that("absorption uses total set abundance and numerical tuple ties", {
  # Lexical string ordering would put 1,10 ahead of 1,2; numerical order must not.
  X <- matrix(0, 5, 12)
  X[1:2, c(1, 10)] <- rep(c(9, 8), each = 2)
  X[3:4, c(1, 2)] <- rep(c(9, 8), each = 2)
  X[5, c(11, 12)] <- c(10, 9)
  zero <- linf.udcsts(X, n0 = 2, low.freq.policy = "absorb")
  expect_identical(zero$retained.ids, c("1,2", "1,10"))
  expect_equal(zero$assignment[5], "1,2")
  X[5, c(1, 2, 10)] <- c(1, 2, 3)
  summed <- linf.udcsts(X, n0 = 2, low.freq.policy = "absorb")
  expect_equal(summed$assignment[5], "1,10")
  # Both candidate sums overflow, but their true ordering differs.
  huge <- matrix(0, 5, 4)
  huge[1:2, c(1, 2)] <- rep(c(9e307, 8e307), each = 2)
  huge[3:4, c(1, 3)] <- rep(c(9e307, 8e307), each = 2)
  huge[5, ] <- c(1e308, 9e307, 1.1e308, 1.2e308)
  expect_equal(linf.udcsts(huge, n0 = 2, low.freq.policy = "absorb")$assignment[5],
               "1,3")
})

test_that("sparse and dense policies agree including zero ranks", {
  X <- rbind(udcst.example(), c(A = 0, B = 0, C = 7))
  rownames(X) <- paste0("sample", seq_len(nrow(X)))
  for (level in 1:3) for (policy in c("pure", "filter", "absorb")) {
    dense <- linf.udcsts(X, level, n0 = 2, low.freq.policy = policy)
    sparse <- linf.udcsts(Matrix::Matrix(X, sparse = TRUE), level,
                          n0 = 2, low.freq.policy = policy)
    expect_identical(sparse$matrix.backend, "sparse")
    sparse$matrix.backend <- "dense"
    expect_identical(sparse, dense)
  }
  expect_equal(linf.udcsts(as.data.frame(X), n0 = 1)$feature.sets,
               linf.udcsts(X, n0 = 1)$feature.sets)
})

test_that("set identity is independent of display-label separators", {
  X <- udcst.example()
  labels <- c("A + B", "1,2", "RARE_DOMINANT")
  fit <- linf.udcsts(X, n0 = 1, feature.ids = c("x", "y", "z"),
                    feature.labels = labels)
  expect_equal(fit$pure.id, c(rep("1,2", 4), "2,3"))
  expect_identical(fit$feature.labels, labels)
  expect_identical(fit$feature.ids, c("x", "y", "z"))
  expect_identical(fit$input.dimnames, dimnames(X))
  expect_error(linf.udcsts(X, feature.ids = "x"), "length ncol")
})

test_that("one row, one feature, and full feature sets retain matrix shape", {
  one <- linf.udcsts(matrix(4, 1, 1), level = 1, n0 = 1)
  expect_identical(dim(one$feature.sets), c(1L, 1L))
  expect_equal(one$pure.id, "1")
  one.row <- linf.udcsts(matrix(c(1, 3, 2), 1), n0 = 1)
  expect_equal(one.row$feature.sets, matrix(c(2L, 3L), 1))
  full <- linf.udcsts(udcst.example(), level = 3, n0 = 5)
  expect_true(all(full$pure.id == "1,2,3"))
  expect_equal(full$frequency.table$freq, 5L)
})

test_that("invalid inputs produce specific errors", {
  X <- udcst.example()
  for (value in list(NA_real_, Inf, 1+1i, 0, 1.5, numeric(), c(1, 2), "2")) {
    expect_error(linf.udcsts(X, level = value), "level must")
    expect_error(linf.udcsts(X, n0 = value), "n0 must")
  }
  expect_error(linf.udcsts(X, level = 4), "level must")
  expect_error(linf.udcsts(matrix(numeric(), 0, 3)), "zero rows or columns")
  expect_error(linf.udcsts(matrix(-1, 1, 2)), "negative")
  expect_error(linf.udcsts(matrix(NA_real_, 1, 2)), "non-finite")
  expect_error(linf.udcsts(matrix(Inf, 1, 2)), "non-finite")
  expect_error(linf.udcsts(matrix("a", 1, 2)), "numeric matrix")
  expect_error(linf.udcsts(1:3), "numeric matrix")
  expect_error(linf.udcsts(X, low.freq.policy = "other"), "arg")
})
