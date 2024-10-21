test_that("print works for NE glycan graphs", {
  x <- new_ne_glycan_graph(n_glycan_core_ne())
  expect_snapshot(print(x))
})


test_that("print works for DN glycan graphs", {
  x <- new_dn_glycan_graph(n_glycan_core_dn())
  expect_snapshot(print(x))
})


test_that("verbose print works for NE glycan graphs", {
  skip_on_old_win()
  x <- new_ne_glycan_graph(n_glycan_core_ne())
  expect_snapshot(print(x, verbose = TRUE))
})


test_that("verbose print works for DN glycan graphs", {
  skip_on_old_win()
  x <- new_dn_glycan_graph(n_glycan_core_dn())
  expect_snapshot(print(x, verbose = TRUE))
})


test_that("verbose print works for NE glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_ne_glycan_graph(n_glycan_core_ne(linkage = FALSE))
  expect_snapshot(print(x, verbose = TRUE))
})


test_that("verbose print works for DN glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_dn_glycan_graph(n_glycan_core_dn(linkage = FALSE))
  expect_snapshot(print(x, verbose = TRUE))
})
