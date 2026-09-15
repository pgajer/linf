test_that("landmarks preserve full IDs instead of matching a suffix feature", {
  M <- rbind(s1 = c(1, 0, .2), s2 = c(.9, .1, .3), s3 = c(.7, .2, .1))
  ids <- c("g__A", "A", "g__B")
  labels <- c("Taxon __ A", "Unrelated suffix", "Taxon __ B")
  for (backend in c("dense", "sparse")) {
    fit <- linf.csts(M, feature.ids = ids, feature.labels = labels,
                     n0 = 2, backend = backend, return.landmarks = TRUE)
    out <- linf.landmarks(M, fit)
    expect_true(out$lineages$landmarks.computable)
    expect_identical(out$landmarks$target.feature.id, rep("g__A", 4))
    expect_equal(out$landmarks$point.name, c("s1", "s3", "s2", "s2"))
    expect_equal(out$landmarks$observed.value, c(1, .7, .9, .9))
    expect_equal(out$landmarks$target.value, c(1, .7, mean(M[, 1]), .9))
    expect_equal(fit$landmarks$landmarks$target.feature.id, rep("g__A", 2))
  }
})

test_that("separator IDs preserve selections through depths and policy views", {
  M <- rbind(
    a1 = c(1, .00, .8, .00, .2, .00),
    a2 = c(.95, .01, .7, .02, .3, .03),
    a3 = c(.9, .03, .2, .02, .7, .01),
    a4 = c(.85, .02, .1, .01, .8, .03),
    b1 = c(.2, .01, 1, .01, .1, .01),
    b2 = c(.3, .01, .95, .01, .1, .01),
    rare = c(.1, .02, .2, .03, 1, .01),
    zero = c(0, 0, 0, 0, 0, 0)
  )
  plain.ids <- c("alpha", "a", "beta", "b", "gamma", "c")
  for (sep in c("__", "::")) {
    ids <- c(paste0("g", sep, "a"), "a", paste0("g", sep, "b"), "b",
             paste0("g", sep, "c"), "c")
    labels <- paste0("Display", sep, seq_along(ids))
    for (backend in c("dense", "sparse")) {
      for (policy in c("pure", "absorb")) {
        plain <- linf.csts(M, feature.ids = plain.ids, feature.labels = labels,
                           n0 = 2, low.freq.policy = policy, backend = backend)
        fit <- linf.csts(M, feature.ids = ids, feature.labels = labels,
                         n0 = 2, low.freq.policy = policy, backend = backend)
        for (depth in 1:3) {
          if (depth > 1) {
            plain <- refine.linf.csts(M, plain, n0 = 2, sep = sep,
                                      low.freq.policy = policy, verbose = FALSE)
            fit <- refine.linf.csts(M, fit, n0 = 2, sep = sep,
                                    low.freq.policy = policy, verbose = FALSE)
          }
          for (view in c("active", "pure", "absorb")) {
            expected <- linf.landmarks(M, plain, depth = depth, view = view)
            actual <- linf.landmarks(M, fit, depth = depth, view = view)
            expect_equal(actual$lineages$landmarks.computable,
                         expected$lineages$landmarks.computable)
            expect_equal(actual$lineages$target.feature.id,
                         ids[match(expected$lineages$target.feature.id, plain.ids)])
            cols <- c("landmark.type", "point.index", "point.name",
                      "target.feature.label", "observed.value", "target.value",
                      "abs.deviation")
            expect_equal(actual$landmarks[, cols], expected$landmarks[, cols])
          }
        }
        # At depth 2, the unrefined beta parent is terminal, not a suffix b.
        terminal <- linf.landmarks(M, fit, depth = 2, view = "pure")$lineages
        expect_equal(terminal$target.feature.id[terminal$lineage.id == ids[3]], ids[3])
      }
    }
  }
})

test_that("rare labels containing separators are skipped at every depth", {
  M <- rbind(a1 = c(10, 5, 1), a2 = c(9, 4, 1), a3 = c(8, 1, 5),
             rare = c(1, 1, 10))
  ids <- c("g__A", "g__B", "g__C")
  for (backend in c("dense", "sparse")) {
    fit <- linf.csts(M, feature.ids = ids, n0 = 2,
                     rare.label = "RARE__BUCKET", backend = backend)
    fit <- refine.linf.csts(M, fit, lineages.to.refine = ids[1], n0 = 2,
                            rare.label = "RARE__BUCKET", verbose = FALSE)
    out <- linf.landmarks(M, fit, view = "pure")
    expect_setequal(out$lineages$lineage.id[out$lineages$is.rare],
                    c("RARE__BUCKET", "g__A__RARE__BUCKET"))
    expect_false(any(out$lineages$landmarks.computable[out$lineages$is.rare]))
    expect_identical(unique(out$landmarks$target.feature.id), "g__B")
  }
})
