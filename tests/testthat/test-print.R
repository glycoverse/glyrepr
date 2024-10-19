test_that("print works for NE glycan graphs", {
  x <- new_ne_glycan_graph(n_glycan_core_ne())
  expect_snapshot(print(x))
})


test_that("print works for DN glycan graphs", {
  x <- new_dn_glycan_graph(n_glycan_core_dn())
  expect_snapshot(print(x))
})
