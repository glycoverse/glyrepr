test_that("print works for glycan graphs", {
  skip_on_old_win()
  x <- new_glycan_graph(n_glycan_core())
  expect_snapshot(print(x))
})


test_that("print works for glycan graphs with verbose = FALSE", {
  x <- new_glycan_graph(n_glycan_core())
  expect_snapshot(print(x, verbose = FALSE))
})


test_that("print works for glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_glycan_graph(n_glycan_core(linkage = FALSE))
  expect_snapshot(print(x))
})


test_that("print works for one-mono glycan graphs", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for graphs with substituent", {
  skip_on_old_win()
  glycan <- o_glycan_core_1()
  igraph::V(glycan)$sub <- c("6S", "")
  expect_snapshot(print(glycan))
})


test_that("print works for one-node graphs with substituent", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- "6S"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for graphs with alditol", {
  skip_on_old_win()
  glycan <- o_glycan_core_1()
  glycan$alditol <- TRUE
  expect_snapshot(print(glycan))
})


test_that("print works for one-node graphs with alditol", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- TRUE
  glycan <- new_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for graph with alditol and substituent on root", {
  skip_on_old_win()
  glycan <- o_glycan_core_1()
  igraph::V(glycan)$sub <- c("6S", "")
  glycan$alditol <- TRUE
  expect_snapshot(print(glycan))
})


test_that("print works for one-node graph with alditol and substituent on root", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- "6S"
  graph$anomer <- "a1"
  graph$alditol <- TRUE
  glycan <- new_glycan_graph(graph)
  expect_snapshot(print(glycan))
})
