test_that("get the composition for a glycan graph", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_glycan_graph(graph)

  comp <- get_composition(glycan)

  expect_equal(sort(comp), sort(c(Glc = 2L, Gal = 1L)))
})
