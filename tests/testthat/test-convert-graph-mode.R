example_ne_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("N", "N", "H")
  igraph::E(graph)$linkage <- c("b1-4", "a1-3")
  new_ne_glycan_graph(graph)
}


example_dn_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 1-+4, 4-+5)
  igraph::V(graph)$type <- c("mono", "linkage", "mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("N", NA, "H", NA, "H")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA, "a1-3", NA)
  new_dn_glycan_graph(graph)
}


test_that("converting NE to DN graph works on simple example", {
  glycan <- example_ne_glycan_graph()
  dn_graph <- convert_ne_to_dn(glycan)
  expect_snapshot(print(dn_graph, verbose = TRUE))
})


test_that("converting NE to DN graph works on complex example", {
  glycan <- n_glycan_core("ne")
  dn_graph <- convert_ne_to_dn(glycan)
  expect_snapshot(print(dn_graph, verbose = TRUE))
})


test_that("converting DN graph fails for DN graphs", {
  glycan <- example_dn_glycan_graph()
  expect_error(convert_ne_to_dn(glycan), "Input must be an NE glycan graph.")
})


test_that("converting one-node graph to DN graph works", {
  graph <- igraph::make_graph(~ 1)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$type <- "mono"
  glycan <- new_ne_glycan_graph(graph)

  dn_graph <- convert_ne_to_dn(glycan)

  expect_snapshot(print(dn_graph, verbose = TRUE))
})
