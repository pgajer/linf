test_that("bundled data provenance strings identify the intended sources", {
  data("valencia2k", package = "linf")
  expect_match(
    valencia2k$source,
    "doi:10[.]1186/s40168-020-00934-6",
    fixed = FALSE
  )
  expect_false(grepl("10[.]1128/mSystems[.]00149-20", valencia2k$source))

  data("agp_gut", package = "linf")
  expect_match(agp_gut$source, "Phenotypes do not influence selection")
  expect_match(agp_gut$source, "selection reasons")
  expect_match(agp_gut$source, "DATA_PROVENANCE[.]md")
  expect_setequal(
    unique(agp_gut$meta$selection_reason),
    c("all_target_dcst", "seeded_random_background")
  )
  expect_equal(
    unname(table(agp_gut$meta$selection_reason)["all_target_dcst"]),
    287L
  )
  expect_equal(
    unname(table(agp_gut$meta$selection_reason)["seeded_random_background"]),
    479L
  )
})

test_that("installed data provenance notice is present", {
  provenance <- system.file("DATA_PROVENANCE.md", package = "linf")
  expect_true(nzchar(provenance))
  expect_true(file.exists(provenance))
})
