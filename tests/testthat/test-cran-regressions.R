test_that("default transfer scores display lineages after ID alignment", {
  X <- rbind(s1 = c(asv1 = 10, asv2 = 2, asv3 = 1),
             s2 = c(9, 3, 1), s3 = c(1, 10, 2), s4 = c(1, 9, 3))
  for (backend in c("dense", "sparse")) {
    M <- normalize.linf(X, backend = backend)
    fit <- linf.csts(M, feature.labels = c("A", "B", "C"),
                     n0 = 2, low.freq.policy = "absorb")
    fit <- refine.linf.csts(M, fit, n0 = 1, low.freq.policy = "absorb",
                            verbose = FALSE)
    for (view in c("absorb", "pure", "active")) {
      tr <- transfer.dcsts(X[, 3:1], fit, view = view, backend = backend)
      expect_equal(tr$assignment[, 1], fit$lineage.labels[[1]])
      expect_equal(tr$assignment[, 2], fit$lineage.labels[[2]])
      expect_equal(tr$assignment, transfer.dcsts(X, fit, view = view,
        match.by = "feature.labels", feature.labels = c("A", "B", "C"),
        backend = backend)$assignment)
    }
  }
})

test_that("unnamed matrices use the fitted synthetic feature IDs", {
  X <- matrix(c(10, 2, 1, 9, 3, 1, 1, 10, 2, 1, 9, 3), ncol = 3, byrow = TRUE)
  for (backend in c("dense", "sparse")) {
    fit <- linf.csts(X, n0 = 2, low.freq.policy = "absorb", backend = backend)
    fit <- refine.linf.csts(X, fit, n0 = 1, low.freq.policy = "absorb", verbose = FALSE)
    expect_equal(transfer.dcsts(X, fit, backend = backend)$assignment[, 1], fit$lineage.labels[[1]])
    expect_equal(transfer.dcsts(X, fit, backend = backend)$assignment[, 2], fit$lineage.labels[[2]])
  }
})

test_that("taxonomy separators do not split feature identities", {
  X <- rbind(s1 = c(10, 3, 1), s2 = c(9, 3, 2), s3 = c(1, 10, 3), s4 = c(2, 9, 3))
  colnames(X) <- c("g__A;s__a", "g__B;s__b", "g__C;s__c")
  for (backend in c("dense", "sparse")) {
    fit <- linf.csts(X, n0 = 2, low.freq.policy = "absorb", backend = backend)
    for (d in 2:3) {
      fit <- refine.linf.csts(X, fit, n0 = 1, low.freq.policy = "absorb", verbose = FALSE)
      expect_equal(transfer.dcsts(X, fit, backend = backend)$assignment[, d], fit$lineage.labels[[d]])
    }
    expect_equal(unname(fit$lineage.ids[[2]][1]), "g__A;s__a__g__B;s__b")
    expect_equal(unname(fit$lineage.ids[[3]][1]), "g__A;s__a__g__B;s__b__g__C;s__c")
  }
})

test_that("transfer never invents an unretained state", {
  X <- rbind(c(A = 10, B = 1, C = 0), c(9, 1, 0), c(1, 10, 0))
  for (backend in c("dense", "sparse")) {
    fit <- linf.csts(X, n0 = 2, low.freq.policy = "absorb", backend = backend)
    query <- rbind(unretained = c(A = 0, B = 0, C = 5), zero = c(0, 0, 0))
    for (carry in c(TRUE, FALSE)) {
      expect_true(all(is.na(transfer.dcsts(query, fit, backend = backend,
        carry.forward.terminal.depths = carry)$assignment)))
      expect_true(all(is.na(transfer.dcsts(query[, "C", drop = FALSE], fit,
        backend = backend, carry.forward.terminal.depths = carry)$assignment)))
      empty <- linf.csts(X, n0 = 4, low.freq.policy = "absorb", backend = backend)
      expect_true(all(is.na(transfer.dcsts(query, empty, backend = backend,
        carry.forward.terminal.depths = carry)$assignment)))
    }
  }
})

test_that("absorption ties and refinement agree across backends", {
  X <- rbind(a1 = c(A = 10, B = 1, C = 0), a2 = c(9, 1, 0),
             b1 = c(1, 10, 0), b2 = c(1, 9, 0), b3 = c(1, 8, 0),
             tied = c(2, 2, 10), zero = c(0, 0, 0))
  for (method in c("first", "random")) {
    fits <- lapply(c("dense", "sparse"), function(backend) {
      set.seed(123)
      linf.csts(X, n0 = 2, low.freq.policy = "absorb",
                tie.method = method, backend = backend)
    })
    expect_equal(fits[[1]]$lineage.labels, fits[[2]]$lineage.labels)
    if (method == "first") expect_equal(unname(fits[[1]]$lineage.label["tied"]), "A")
    expect_equal(unname(fits[[1]]$lineage.label["zero"]), "B")
    refined <- lapply(seq_along(fits), function(i) {
      set.seed(123)
      refine.linf.csts(X, fits[[i]], n0 = 1, low.freq.policy = "absorb",
                       verbose = FALSE)
    })
    expect_equal(refined[[1]]$lineage.labels, refined[[2]]$lineage.labels)
  }
  for (backend in c("dense", "sparse")) {
    expect_error(linf.csts(X, n0 = 2, low.freq.policy = "absorb",
                          tie.method = "error", backend = backend), "tie during reassignment")
  }
})

test_that("feature metadata lengths are validated for both backends", {
  X <- matrix(1:3, nrow = 1)
  for (backend in c("dense", "sparse")) {
    for (fun in list(linf.dominant.features, linf.csts)) {
      expect_error(fun(X, feature.ids = "only_one", backend = backend), "ncol")
      expect_error(fun(X, feature.labels = "only_one", backend = backend), "length")
    }
  }
})

test_that("normalization leaves below-tolerance rows unchanged", {
  X <- rbind(tiny = c(A = .01, B = .005), zero = c(0, 0))
  for (backend in c("dense", "sparse")) {
    expect_equal(as.matrix(normalize.linf(X, tol = .1, backend = backend)), X)
    expect_equal(unname(linf.dominant.features(X, backend = backend)$label), c("A", NA))
    for (tol in list(NA_real_, NaN, Inf, -1, numeric())) {
      expect_error(normalize.linf(X, tol = tol, backend = backend), "tol must")
    }
  }
})
