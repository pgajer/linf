test_that("normalize.linf max-per-row", {
  S <- rbind(c(1, 2, 4), c(0, 0, 0))
  Z <- normalize.linf(S)
  expect_equal(apply(Z, 1, max), c(1, 0))
})

test_that("linf.dominant.features ties, zeros, labels, levels", {
  S <- rbind(c(0.7, 0.3),
             c(0, 0),
             c(1, 1))   # tie -> first
  colnames(S) <- c("A", "B")

  out <- linf.dominant.features(S)

  # index behavior
  expect_equal(out$index, c(1L, NA_integer_, 1L))

  # label behavior (mapped from levels)
  expect_equal(out$label, c("A", NA_character_, "A"))

  # levels must be full column-ordered set, 1–1 with columns
  expect_equal(out$levels, c("A", "B"))

  # observed.levels is a subset of levels, in column order, no NA
  expect_equal(out$observed.levels, "A")
})

test_that("linf.dominant.features unnamed columns synthesize unique levels", {
  T <- matrix(c(0,2,  3,1,  0,0), nrow = 3, byrow = TRUE)
  # no colnames -> V1, V2
  out <- linf.dominant.features(T)

  expect_equal(out$levels, c("V1", "V2"))     # synthesized
  expect_equal(out$index, c(2L, 1L, NA_integer_))
  expect_equal(out$label, c("V2", "V1", NA_character_))
  expect_equal(out$observed.levels, c("V1", "V2"))  # both appeared
})

test_that("linf.dominant.features duplicate column names get disambiguated", {
  T <- matrix(c(0, 2,
                3, 1,
                0, 0), nrow = 3, byrow = TRUE)
  colnames(T) <- c("X", "X") # duplicates

  out <- linf.dominant.features(T)
  # expected: make.unique -> "X", "X_1"
  expect_equal(out$levels, c("X", "X_1"))
  expect_equal(out$label[1:2], c("X_1", "X"))
})

test_that("linf.csts retains supported states and reassigns low-support ones", {
  # 6 rows dominated by col1, 2 by col2, 2 by col3
  A <- matrix(c(5, 1, 0), nrow = 6, ncol = 3, byrow = TRUE)
  B <- matrix(c(1, 5, 0), nrow = 2, ncol = 3, byrow = TRUE)
  C <- matrix(c(1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
  S <- rbind(A, B, C)
  S <- sweep(S, 1, rowSums(S), "/")
  colnames(S) <- c("Dom1", "Dom2", "Dom3")

  out <- linf.csts(S, n0 = 5, low.freq.policy = "absorb")

  # only col1 kept
  expect_identical(out$retained.feature.indices, 1L)
  expect_identical(out$retained.feature.labels, "Dom1")

  # supported state unchanged
  expect_true(all(out$depth1.feature.index[1:6] == 1L))
  expect_true(all(out$lineage.label[1:6] == "Dom1"))

  # low-support samples reassigned to the retained state
  expect_true(all(out$depth1.feature.index[7:10] == 1L))
  expect_true(all(out$lineage.label[7:10] == "Dom1"))
})

test_that("linf.csts returns all-NA when no states meet the support threshold", {
  S <- rbind(
    c(1, 0, 0),  # 1 sample in col1
    c(0, 1, 0),  # 1 sample in col2
    c(0, 0, 1)   # 1 sample in col3
  )
  colnames(S) <- c("A","B","C")
  out <- linf.csts(S, n0 = 2)

  expect_true(all(is.na(out$depth1.feature.index)))
  expect_true(all(out$lineage.label == out$rare.label))
  expect_length(out$retained.feature.indices, 0)
  expect_length(out$retained.feature.labels, 0)

  # size.table reports raw sizes by label
  expect_identical(as.integer(unclass(out$size.table)), c(1L,1L,1L))
  expect_identical(names(out$size.table), c("A","B","C"))
})
