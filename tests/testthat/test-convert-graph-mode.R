simple_example_ne_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("N", "N", "H")
  igraph::E(graph)$linkage <- c("b1-4", "a1-3")
  new_ne_glycan_graph(graph)
}


simple_example_dn_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 1-+4, 4-+5)
  igraph::V(graph)$type <- c("mono", "linkage", "mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("N", NA, "H", NA, "H")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA, "a1-3", NA)
  new_dn_glycan_graph(graph)
}


# This is actually an N-glycan core.
# We don't use `n_glycan_core()` here because it actually depends on
# `convert_ne_to_dn()` to create a DN graph.
# To prevent circular dependency
# (it will not raise an error but the tests will always pass),
# we create the graph manually here.
complex_example_ne_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4", "a1-3", "a1-6")
  new_ne_glycan_graph(graph)
}


complex_example_dn_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 4-+5, 5-+6, 6-+7, 5-+8, 8-+9)
  igraph::V(graph)$type <- c("mono", rep(c("linkage", "mono"), times = 4))
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc", NA, "Man", NA, "Man", NA, "Man")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA, "b1-4", NA, "a1-3", NA, "a1-6", NA)
  new_dn_glycan_graph(graph)
}


test_that("converting NE to DN graph works on simple example", {
  skip_on_old_win()
  glycan <- simple_example_ne_glycan_graph()
  dn_graph <- convert_ne_to_dn(glycan)
  expect_snapshot(print(dn_graph, verbose = TRUE))
})


test_that("converting NE to DN graph works on complex example", {
  skip_on_old_win()
  glycan <- complex_example_ne_glycan_graph()
  dn_graph <- convert_ne_to_dn(glycan)
  expect_snapshot(print(dn_graph, verbose = TRUE))
})



test_that("converting to DN graph fails for DN graphs", {
  glycan <- simple_example_dn_glycan_graph()
  expect_error(convert_ne_to_dn(glycan), "Input must be an NE glycan graph.")
})


test_that("converting one-node graph to DN graph works", {
  skip_on_old_win()
  graph <- igraph::make_graph(~ 1)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$type <- "mono"
  glycan <- new_ne_glycan_graph(graph)

  dn_graph <- convert_ne_to_dn(glycan)

  expect_snapshot(print(dn_graph, verbose = TRUE))
})


test_that("converting DN to EN graph works on simple example", {
  skip_on_old_win()
  glycan <- simple_example_dn_glycan_graph()
  ne_graph <- convert_dn_to_ne(glycan)
  expect_snapshot(print(ne_graph, verbose = TRUE))
})


test_that("converting DN to EN graph works on complex example", {
  skip_on_old_win()
  glycan <- complex_example_dn_glycan_graph()
  ne_graph <- convert_dn_to_ne(glycan)
  expect_snapshot(print(ne_graph, verbose = TRUE))
})


test_that("converting to NE graph fails for NE graphs", {
  glycan <- simple_example_ne_glycan_graph()
  expect_error(convert_dn_to_ne(glycan), "Input must be a DN glycan graph.")
})


test_that("converting one-node graph to NE graph works", {
  skip_on_old_win()
  graph <- igraph::make_graph(~ 1)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$linkage <- NA_character_
  glycan <- new_dn_glycan_graph(graph)

  ne_graph <- convert_dn_to_ne(glycan)

  expect_snapshot(print(ne_graph, verbose = TRUE))
})


test_that("ensureing NE glycan as NE", {
  glycan <- simple_example_ne_glycan_graph()
  glycan <- ensure_graph_mode(glycan, "ne")
  expect_true(is_ne_glycan(glycan))
})


test_that("ensureing NE glycan as DN", {
  glycan <- simple_example_ne_glycan_graph()
  glycan <- ensure_graph_mode(glycan, "dn")
  expect_true(is_dn_glycan(glycan))
})


test_that("ensureing DN glycan as DN", {
  glycan <- simple_example_dn_glycan_graph()
  glycan <- ensure_graph_mode(glycan, "dn")
  expect_true(is_dn_glycan(glycan))
})


test_that("ensureing DN glycan as NE", {
  glycan <- simple_example_dn_glycan_graph()
  glycan <- ensure_graph_mode(glycan, "ne")
  expect_true(is_ne_glycan(glycan))
})
