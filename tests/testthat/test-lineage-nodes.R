collision_matrix <- function() rbind(
  a1 = c(A = 10, B = 5, A__B = 1), a2 = c(A = 9, B = 4, A__B = 1),
  ab1 = c(A = 1, B = 0, A__B = 10), ab2 = c(A = 1, B = 0, A__B = 9)
)

test_that("explicit nodes separate a literal feature from a refined path", {
  X <- collision_matrix()
  for (backend in c("dense", "sparse")) for (policy in c("pure", "absorb")) {
    M <- if (backend == "sparse") Matrix::Matrix(X, sparse = TRUE) else X
    f <- linf.csts(M, n0 = 2, low.freq.policy = policy)
    f <- refine.linf.csts(M, f, lineages.to.refine = "A", n0 = 2,
                          low.freq.policy = policy, verbose = FALSE)
    expect_length(unique(f$lineage.id), 2L)
    expect_length(unique(f$lineage.label), 2L)
    for (view in c("active", "pure", "absorb")) {
      lm <- linf.landmarks(M, f, view = view)
      expect_equal(lm$lineages$lineage.size, c(2L, 2L))
      expect_equal(lm$lineages$target.feature.id, c("B", "A__B"))
      expect_equal(unname(transfer.dcsts(M, f, view = view)$assignment[, 2]),
                   unname(f$lineage.label))
    }
    # Further refinement removes A and B from their path, but only the literal
    # A__B feature from the terminal path.
    f3 <- refine.linf.csts(M, f, lineages.to.refine = unique(f$lineage.id),
                           n0 = 2, low.freq.policy = policy, verbose = FALSE)
    expect_equal(linf.landmarks(M, f3)$lineages$target.feature.id, c("A__B", "A"))
    expect_equal(unname(transfer.dcsts(M, f3)$assignment[, 3]), unname(f3$lineage.label))
  }
})

test_that("colliding display labels and synthetic names retain distinct identities", {
  X <- collision_matrix()
  colnames(X) <- c("id1", "id2", "id3")
  f <- linf.csts(X, feature.labels = c("A", "B", "A::B"), n0 = 2)
  f <- refine.linf.csts(X, f, lineages.to.refine = "id1", sep = "::", n0 = 2, verbose = FALSE)
  expect_length(unique(f$lineage.label), 2L)
  expect_equal(unname(transfer.dcsts(X, f)$assignment[, 2]), unname(f$lineage.label))
  R <- rbind(c(10, 1), c(10, 2), c(1, 10))
  for (rare in c("RARE_DOMINANT", "rare::state")) {
    colnames(R) <- c(rare, "B")
    f <- linf.csts(R, n0 = 2, rare.label = rare)
    lm <- linf.landmarks(R, f)
    expect_equal(sort(lm$lineages$lineage.size), c(1L, 2L))
    expect_equal(sum(lm$lineages$is.rare), 1L)
    expect_equal(unique(lm$landmarks$target.feature.id), rare)
    expect_equal(summary(f)$rare.samples, 1L)
    expect_equal(unname(transfer.dcsts(R[1:2, , drop = FALSE], f, view = "pure")$assignment[, 1]),
                 unname(f$lineage.label[1:2]))
    real.id <- lm$lineages$lineage.id[!lm$lineages$is.rare]
    expect_no_error(refine.linf.csts(R, f, lineages.to.refine = real.id,
                                    rare.label = rare, n0 = 1, verbose = FALSE))
  }
})

test_that("node assignments survive row permutation and RDS serialization", {
  X <- collision_matrix()
  fit <- function(M) refine.linf.csts(M, linf.csts(M, n0 = 2),
                                     lineages.to.refine = "A", n0 = 2, verbose = FALSE)
  f <- fit(X)
  p <- c(4, 2, 3, 1)
  g <- fit(X[p, ])
  expect_equal(g$lineage.id[order(p)], f$lineage.id)
  tmp <- tempfile(fileext = ".rds")
  on.exit(unlink(tmp))
  saveRDS(f, tmp)
  expect_equal(linf.landmarks(X, readRDS(tmp)), linf.landmarks(X, f))
  # Nodes are based on feature paths, not a split of displayed IDs.
  expect_identical(f$nodes[[2]]$path[[1]], c(1L, 2L))
  expect_true(f$nodes[[2]]$terminal[f$nodes[[2]]$node.id == "3"])
})

test_that("unambiguous legacy fits upgrade and ambiguous legacy fits fail safely", {
  legacy <- function(f) { f[grep("^(nodes|node.ids)", names(f))] <- NULL; f }
  X <- collision_matrix(); colnames(X)[3] <- "AB"
  f <- refine.linf.csts(X, linf.csts(X, n0 = 2), lineages.to.refine = "A", n0 = 2, verbose = FALSE)
  expect_equal(linf.landmarks(X, legacy(f)), linf.landmarks(X, f))
  expect_equal(transfer.dcsts(X, legacy(f)), transfer.dcsts(X, f))
  # Reconstruct the historically ambiguous strings, retaining their levels.
  old <- legacy(f)
  old$feature.ids[3] <- old$feature.labels[3] <- "A__B"
  for (field in grep("^lineage", names(old), value = TRUE)) {
    if (is.list(old[[field]])) old[[field]] <- lapply(old[[field]], function(x) replace(x, x == "AB", "A__B"))
    else old[[field]][old[[field]] == "AB"] <- "A__B"
  }
  for (consume in list(function() linf.landmarks(X, old), function() transfer.dcsts(X, old),
                       function() refine.linf.csts(X, old, verbose = FALSE), function() summary(old))) {
    expect_error(consume(), "ambiguous legacy.*rebuild")
  }
})


test_that("synthetic parents remain terminal when the absorbed view is refined", {
  X <- rbind(a1 = c(A = 10, B = 5, C = 1),
             a2 = c(A = 9, B = 4, C = 1),
             rare = c(A = 8, B = 1, C = 10))
  fit <- function(M) refine.linf.csts(M, linf.csts(M, n0 = 2, low.freq.policy = "absorb"),
                                     lineages.to.refine = "A", n0 = 2,
                                     low.freq.policy = "absorb", verbose = FALSE)
  f <- fit(X)
  expect_equal(unname(f$lineage.id.pure), c("A__B", "A__B", "RARE_DOMINANT"))
  expect_equal(unname(f$lineage.id.absorb), rep("A__B", 3))
  expect_equal(summary(dcst.view(f, "pure"))$rare.samples, c(1L, 1L))
  p <- c(3, 2, 1)
  expect_equal(fit(X[p, ])$lineage.id.pure[order(p)], f$lineage.id.pure)
})
