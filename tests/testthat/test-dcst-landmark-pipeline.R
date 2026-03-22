test_that("linf.dcst.landmark.pipeline returns the expected object bundle", {
  X <- rbind(
    s1 = c(10, 8, 1),
    s2 = c(9, 7, 2),
    s3 = c(8, 2, 7),
    s4 = c(7, 1, 8),
    s5 = c(1, 10, 2),
    s6 = c(2, 9, 1)
  )
  ids <- c("asv_1", "asv_2", "asv_3")
  labels <- c("L. iners 1", "Gard. vaginalis 2", "BVAB1 3")

  out <- linf.dcst.landmark.pipeline(
    X,
    feature.ids = ids,
    feature.labels = labels,
    n0.depth1 = 2,
    n0.depth2 = 2,
    refinement.factor = 2,
    low.freq.policy = "rare",
    landmark.view = "absorb",
    verbose = FALSE
  )

  expect_named(
    out,
    c("linf.rel", "dcst.depth1", "dcst.depth2", "landmarks.depth1", "landmarks.depth2", "params")
  )
  expect_equal(dim(out$linf.rel), dim(X))
  expect_s3_class(out$dcst.depth1, "linf.csts")
  expect_s3_class(out$dcst.depth2, "linf.csts")
  expect_s3_class(out$landmarks.depth1, "linf.landmarks")
  expect_s3_class(out$landmarks.depth2, "linf.landmarks")
  expect_equal(out$dcst.depth1$cst.depth, 1L)
  expect_equal(out$dcst.depth2$cst.depth, 2L)
})

test_that("linf.dcst.landmark.pipeline computes landmark tables for both depths", {
  X <- rbind(
    r1 = c(1.00, 0.80, 0.20),
    r2 = c(0.95, 0.70, 0.30),
    r3 = c(0.90, 0.20, 0.70),
    r4 = c(0.85, 0.10, 0.80),
    r5 = c(0.20, 1.00, 0.10),
    r6 = c(0.30, 0.95, 0.10)
  )
  ids <- c("asv_1", "asv_4", "asv_6")
  labels <- c("L. iners 1", "L. iners 4", "BVAB1 6")

  out <- linf.dcst.landmark.pipeline(
    X,
    feature.ids = ids,
    feature.labels = labels,
    n0.depth1 = 2,
    n0.depth2 = 2,
    refinement.factor = 2,
    low.freq.policy = "rare",
    landmark.view = "absorb",
    depth1.landmark.types = c("endpoint.max", "mean.rep"),
    depth2.landmark.types = c("endpoint.max", "median.rep"),
    verbose = FALSE
  )

  expect_true(all(c("cells", "landmarks") %in% names(out$landmarks.depth1)))
  expect_true(all(c("cells", "landmarks") %in% names(out$landmarks.depth2)))
  expect_gt(nrow(out$landmarks.depth1$cells), 0L)
  expect_gt(nrow(out$landmarks.depth1$landmarks), 0L)
  expect_gt(nrow(out$landmarks.depth2$cells), 0L)
  expect_gt(nrow(out$landmarks.depth2$landmarks), 0L)
  expect_equal(sort(unique(out$landmarks.depth1$landmarks$landmark.type)), c("endpoint.max", "mean.rep"))
  expect_equal(sort(unique(out$landmarks.depth2$landmarks$landmark.type)), c("endpoint.max", "median.rep"))
})

test_that("linf.dcst.landmark.pipeline respects absorb landmark views under rare policy", {
  A <- matrix(c(5, 1, 0), nrow = 6, ncol = 3, byrow = TRUE)
  B <- matrix(c(1, 5, 0), nrow = 2, ncol = 3, byrow = TRUE)
  C <- matrix(c(1, 0, 5), nrow = 2, ncol = 3, byrow = TRUE)
  X <- rbind(A, B, C)
  rownames(X) <- paste0("s", seq_len(nrow(X)))
  colnames(X) <- c("Dom1", "Dom2", "Dom3")

  out <- linf.dcst.landmark.pipeline(
    X,
    n0.depth1 = 5,
    n0.depth2 = 3,
    refinement.factor = 2,
    low.freq.policy = "rare",
    landmark.view = "absorb",
    depth1.landmark.types = "endpoint.max",
    depth2.landmark.types = "endpoint.max",
    verbose = FALSE
  )

  expect_identical(out$dcst.depth1$low.freq.policy, "rare")
  expect_identical(out$landmarks.depth1$view, "absorb")
  expect_identical(out$landmarks.depth2$view, "absorb")
  expect_false(any(out$landmarks.depth1$cells$cell.id == out$dcst.depth1$rare.label))
  expect_false(any(out$landmarks.depth2$cells$cell.id == out$dcst.depth2$rare.label))
})

test_that("linf.dcst.landmark.pipeline carries feature ids and labels into landmark output", {
  X <- rbind(
    s1 = c(0.95, 0.80, 0.05),
    s2 = c(0.90, 0.75, 0.10),
    s3 = c(0.85, 0.20, 0.75),
    s4 = c(0.84, 0.15, 0.80),
    s5 = c(0.10, 0.95, 0.05),
    s6 = c(0.05, 0.94, 0.08)
  )
  ids <- c("asv_1", "asv_4", "asv_5")
  labels <- c("L. iners 1", "L. iners 4", "Mega. lornae")

  out <- linf.dcst.landmark.pipeline(
    X,
    feature.ids = ids,
    feature.labels = labels,
    n0.depth1 = 2,
    n0.depth2 = 2,
    refinement.factor = 2,
    low.freq.policy = "rare",
    landmark.view = "absorb",
    depth1.landmark.types = "endpoint.max",
    depth2.landmark.types = "endpoint.max",
    verbose = FALSE
  )

  expect_true(all(out$dcst.depth1$feature.ids == ids))
  expect_true(all(out$dcst.depth1$feature.labels == labels))
  expect_true(all(out$dcst.depth2$feature.ids == ids))
  expect_true(all(out$dcst.depth2$feature.labels == labels))
  expect_true(all(out$landmarks.depth1$cells$target.feature.id %in% c(ids, NA_character_)))
  expect_true(all(out$landmarks.depth1$cells$target.feature.label %in% c(labels, NA_character_)))
  expect_true(any(grepl("^asv_1(__|$)", out$landmarks.depth2$cells$cell.id)))
  expect_true(any(grepl("^L\\. iners 1(__|$)", out$landmarks.depth2$cells$cell.label)))
})
