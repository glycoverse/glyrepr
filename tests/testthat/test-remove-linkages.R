test_that("works for glycans", {
  glycan <- o_glycan_core_2(linkage = TRUE)
  glycan <- remove_linkages(glycan)
  expect_false(has_linkages(glycan))
})
