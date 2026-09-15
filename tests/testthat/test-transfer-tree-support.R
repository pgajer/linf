test_that("frozen support counts include only observed parent-child pairs", {
  levels <- list(c("p2", "p1", "p1", "p2", "p1", NA, "p3"),
                 c("c3", "c2", "c1", "c3", "c2", "c1", NA))
  nodes <- list(
    data.frame(node.id = c("p1", "p2", "p3"), lineage.label = c("A", "B", "C"),
               feature.index = 1:3, is.rare = FALSE),
    data.frame(node.id = c("c1", "c2", "c3"), lineage.label = c("AA", "AB", "BA"),
               feature.index = c(2L, 3L, 1L), is.rare = FALSE)
  )
  tree <- linf.dcst.transfer.tree(levels, nodes, 2L, "__")
  expect_identical(tree$children[[2]]$p1, c("c2", "c1"))
  expect_identical(tree$children[[2]]$p2, "c3")
  expect_null(tree$children[[2]]$p3)
  expect_equal(unname(tree$support[[2]]$p1[c("c1", "c2")]), c(1, 2))
  expect_equal(unname(tree$support[[2]]$p2["c3"]), 2)
  expect_identical(tree$feature.indices[[2]]$p1, c(3L, 2L))
})
