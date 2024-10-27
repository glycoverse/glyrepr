test_that("works for NE glycans", {
  glycan <- o_glycan_core_2(mode = "ne", linkage = TRUE)
  glycan <- remove_linkages(glycan)
  expect_false(has_linkages(glycan))
})


test_that("works for DN glycans", {
  glycan <- o_glycan_core_2(mode = "dn", linkage = TRUE)
  glycan <- remove_linkages(glycan)
  expect_false(has_linkages(glycan))
})
