test_that("is_glycan works for glycan graphs", {
  x <- new_ne_glycan_graph(igraph::make_empty_graph())
  expect_true(is_glycan(x))
})


test_that("is_glycan not works for others", {
  x <- igraph::make_empty_graph()
  expect_false(is_glycan(x))
})


test_that("is_ne_glycan works", {
  x <- new_ne_glycan_graph(igraph::make_empty_graph())
  expect_true(is_ne_glycan(x))
})


test_that("is_dn_glycan works", {
  x <- new_dn_glycan_graph(igraph::make_empty_graph())
  expect_true(is_dn_glycan(x))
})


test_that("glycan_graph_mode works for NE glycans", {
  glycan <- o_glycan_core_1(mode = "ne")
  expect_equal(glycan_graph_mode(glycan), "ne")
})


test_that("glycan_graph_mode works for DN glycans", {
  glycan <- o_glycan_core_1(mode = "dn")
  expect_equal(glycan_graph_mode(glycan), "dn")
})
