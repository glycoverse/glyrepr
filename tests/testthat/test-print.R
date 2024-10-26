test_that("print works for NE glycan graphs", {
  skip_on_old_win()
  x <- new_ne_glycan_graph(n_glycan_core(mode = "ne"))
  expect_snapshot(print(x))
})


test_that("print works for DN glycan graphs", {
  skip_on_old_win()
  x <- new_dn_glycan_graph(n_glycan_core(mode = "dn"))
  expect_snapshot(print(x))
})


test_that("print works for NE glycan graphs with verbose = FALSE", {
  x <- new_ne_glycan_graph(n_glycan_core(mode = "ne"))
  expect_snapshot(print(x, verbose = FALSE))
})


test_that("print works for DN glycan graphs with verbose = FALSE", {
  x <- new_dn_glycan_graph(n_glycan_core(mode = "dn"))
  expect_snapshot(print(x, verbose = FALSE))
})


test_that("print works for NE glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_ne_glycan_graph(n_glycan_core(mode = "ne", linkage = FALSE))
  expect_snapshot(print(x))
})


test_that("print works for DN glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_dn_glycan_graph(n_glycan_core(mode = "dn", linkage = FALSE))
  expect_snapshot(print(x))
})


test_that("print works for one-mono NE glycan graphs", {
  skip_on_old_win()
  graph <- igraph::make_graph(~ 1)
  igraph::V(graph)$mono <- "N"
  glycan <- new_ne_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for one-mono DN glycan graphs", {
  skip_on_old_win()
  graph <- igraph::make_graph(~ 1)
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$linkage <- NA_character_
  glycan <- new_dn_glycan_graph(graph)
  expect_snapshot(print(glycan))
})
