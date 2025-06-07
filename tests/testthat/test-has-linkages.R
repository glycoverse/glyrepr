test_that("`has_linkages()` works for glycans with intact linkages", {
  glycan <- o_glycan_core_1(linkage = TRUE)
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for glycan with partial linkages", {
  glycan <- o_glycan_core_2(linkage = TRUE)
  glycan <- igraph::set_edge_attr(glycan, "linkage", value = c("??-?", "b1-6"))
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for glycans without linkages", {
  glycan <- o_glycan_core_1(linkage = FALSE)
  expect_false(has_linkages(glycan))
})
