test_that("get_anomer works", {
  x <- n_glycan_core()
  expect_equal(get_anomer(x), "b1")
})

test_that("get_anomer works with a glycan graph", {
  graph <- get_structure_graphs(n_glycan_core())

  expect_identical(get_anomer(graph), "b1")
})
