test_that("normalize.linf max-per-row", {
  S <- rbind(c(1, 2, 4), c(0, 0, 0))
  Z <- normalize.linf(S)
  expect_equal(apply(Z, 1, max), c(1, 0))
})

test_that("linf.cells ties, zeros, labels, levels", {
  S <- rbind(c(0.7, 0.3),
             c(0, 0),
             c(1, 1))   # tie -> first
  colnames(S) <- c("A", "B")

  out <- linf.cells(S)

  # index behavior
  expect_equal(out$index, c(1L, NA_integer_, 1L))

  # label behavior (mapped from levels)
  expect_equal(out$label, c("A", NA_character_, "A"))

  # levels must be full column-ordered set, 1–1 with columns
  expect_equal(out$levels, c("A", "B"))

  # observed.levels is a subset of levels, in column order, no NA
  expect_equal(out$observed.levels, "A")
})

test_that("linf.cells unnamed columns synthesize unique levels", {
  T <- matrix(c(0,2,  3,1,  0,0), nrow = 3, byrow = TRUE)
  # no colnames -> V1, V2
  out <- linf.cells(T)

  expect_equal(out$levels, c("V1", "V2"))     # synthesized
  expect_equal(out$index, c(2L, 1L, NA_integer_))
  expect_equal(out$label, c("V2", "V1", NA_character_))
  expect_equal(out$observed.levels, c("V1", "V2"))  # both appeared
})

test_that("linf.cells duplicate column names get disambiguated", {
  T <- matrix(c(0, 2,
                3, 1,
                0, 0), nrow = 3, byrow = TRUE)
  colnames(T) <- c("X", "X") # duplicates

  out <- linf.cells(T)
  # expected: make.unique -> "X", "X_1"
  expect_equal(out$levels, c("X", "X_1"))
  expect_equal(out$label[1:2], c("X_1", "X"))
})

test_that("linf.csts keeps big cells and reassigns small ones", {
  # 6 rows dominated by col1, 2 by col2, 2 by col3
  A <- matrix(c(5, 1, 0), nrow = 6, ncol = 3, byrow = TRUE)
  B <- matrix(c(1, 5, 0), nrow = 2, ncol = 3, byrow = TRUE)
  C <- matrix(c(1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
  S <- rbind(A, B, C)
  S <- sweep(S, 1, rowSums(S), "/")
  colnames(S) <- c("Dom1", "Dom2", "Dom3")

  out <- linf.csts(S, n0 = 5)

  # only col1 kept
  expect_identical(out$kept.cells.idx, 1L)
  expect_identical(out$kept.cells.lbl, "Dom1")

  # big cell unchanged
  expect_true(all(out$cell.index[1:6] == 1L))
  expect_true(all(out$cell.label[1:6] == "Dom1"))

  # small cells reassigned to kept cell
  expect_true(all(out$cell.index[7:10] == 1L))
  expect_true(all(out$cell.label[7:10] == "Dom1"))
})

test_that("linf.csts returns all-NA when no cells meet threshold", {
  S <- rbind(
    c(1, 0, 0),  # 1 sample in col1
    c(0, 1, 0),  # 1 sample in col2
    c(0, 0, 1)   # 1 sample in col3
  )
  colnames(S) <- c("A","B","C")
  out <- linf.csts(S, n0 = 2)

  expect_true(all(is.na(out$cell.index)))
  expect_true(all(is.na(out$cell.label)))
  expect_length(out$kept.cells.idx, 0)
  expect_length(out$kept.cells.lbl, 0)

  # size.table reports raw sizes by label
  expect_identical(as.integer(unclass(out$size.table)), c(1L,1L,1L))
  expect_identical(names(out$size.table), c("A","B","C"))
})
