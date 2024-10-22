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
