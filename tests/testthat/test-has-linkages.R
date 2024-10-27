test_that("`has_linkages()` works for NE glycans with intact linkages", {
  glycan <- o_glycan_core_1(mode = "ne", linkage = TRUE)
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for NE glycan with partial linakges", {
  glycan <- o_glycan_core_2(mode = "ne", linkage = TRUE)
  glycan <- igraph::set_edge_attr(glycan, "linkage", value = c("??-?", "b1-6"))
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for NE glycans without linkages", {
  glycan <- o_glycan_core_1(mode = "ne", linkage = FALSE)
  expect_false(has_linkages(glycan))
})


test_that("`has_linkages()` works for DN glycans with intact linkages", {
  glycan <- o_glycan_core_1(mode = "dn", linkage = TRUE)
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for DN glycan with partial linakges", {
  glycan <- o_glycan_core_2(mode = "dn", linkage = TRUE)
  glycan <- igraph::set_vertex_attr(glycan, "linkage", value = c(NA, NA, NA, "??-?", "b1-6"))
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for DN glycans without linkages", {
  glycan <- o_glycan_core_1(mode = "dn", linkage = FALSE)
  expect_false(has_linkages(glycan))
})
