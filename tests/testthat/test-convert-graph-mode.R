# !!!!! IMPORTANT !!!!!
# We don't use `n_glycan_core()` or `o_glycan_core_1()` here because it actually
# depends on `convert_ne_to_dn()` to create a DN graph.
# To prevent circular dependency
# (it will not raise an error but the tests will always pass),
# we create the graph manually here.

simple_example_ne_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("N", "N", "H")
  igraph::V(graph)$sub <- c("6S", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "a1-3")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  new_ne_glycan_graph(graph)
}


simple_example_dn_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 1-+4, 4-+5)
  igraph::V(graph)$type <- c("mono", "linkage", "mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("N", NA, "H", NA, "H")
  igraph::V(graph)$sub <- c("6S", NA, "", NA, "")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA, "a1-3", NA)
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  new_dn_glycan_graph(graph)
}


complex_example_ne_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  igraph::V(graph)$sub <- c("6S", "", "", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4", "a1-3", "a1-6")
  graph$anomer <- "?1"
  graph$alditol <- FALSE
  new_ne_glycan_graph(graph)
}


complex_example_dn_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 4-+5, 5-+6, 6-+7, 5-+8, 8-+9)
  igraph::V(graph)$type <- c("mono", rep(c("linkage", "mono"), times = 4))
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc", NA, "Man", NA, "Man", NA, "Man")
  igraph::V(graph)$sub <- c("6S", NA, "", NA, "", NA, "", NA, "")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA, "b1-4", NA, "a1-3", NA, "a1-6", NA)
  graph$anomer <- "?1"
  graph$alditol <- FALSE
  new_dn_glycan_graph(graph)
}


test_that("converting NE to DN graph works on simple example", {
  skip_on_old_win()
  glycan <- simple_example_ne_glycan_graph()
  dn_graph <- convert_graph_mode(glycan, "dn")
  expect_snapshot(print(dn_graph, verbose = TRUE))
  expect_no_error(validate_dn_glycan_graph(dn_graph))
})


test_that("converting NE to DN graph works on complex example", {
  skip_on_old_win()
  glycan <- complex_example_ne_glycan_graph()
  dn_graph <- convert_graph_mode(glycan, "dn")
  expect_snapshot(print(dn_graph, verbose = TRUE))
  expect_no_error(validate_dn_glycan_graph(dn_graph))
})


test_that("converting to DN graph fails for DN graphs", {
  glycan <- simple_example_dn_glycan_graph()
  expect_snapshot(convert_graph_mode(glycan, "dn"), error = TRUE)
})


test_that("converting from DN to DN returns self with `strict` to `FALSE`", {
  glycan <- simple_example_dn_glycan_graph()
  dn_graph <- convert_graph_mode(glycan, "dn", strict = FALSE)
  expect_identical(glycan, dn_graph)
  expect_no_error(validate_dn_glycan_graph(dn_graph))
})


test_that("converting one-node graph to DN graph works", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  igraph::V(graph)$type <- "mono"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_ne_glycan_graph(graph)

  dn_graph <- convert_graph_mode(glycan, "dn")

  expect_snapshot(print(dn_graph, verbose = TRUE))
  expect_no_error(validate_dn_glycan_graph(dn_graph))
})


test_that("converting DN to NE graph works on simple example", {
  skip_on_old_win()
  glycan <- simple_example_dn_glycan_graph()
  ne_graph <- convert_graph_mode(glycan, "ne")
  expect_snapshot(print(ne_graph, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(ne_graph))
})


test_that("converting DN to NE graph works on complex example", {
  skip_on_old_win()
  glycan <- complex_example_dn_glycan_graph()
  ne_graph <- convert_graph_mode(glycan, "ne")
  expect_snapshot(print(ne_graph, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(ne_graph))
})


test_that("converting one-node graph to NE graph works", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$linkage <- NA_character_
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_dn_glycan_graph(graph)

  ne_graph <- convert_graph_mode(glycan, "ne")

  expect_snapshot(print(ne_graph, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(ne_graph))
})
