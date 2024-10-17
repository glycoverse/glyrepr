test_that("glycan graph class", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")

  glycan <- new_glycan_graph(graph)

  expect_s3_class(glycan, c("glycan_graph", "igraph"))
})


test_that("validating undirected graphs", {
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc")
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must be directed")
})


test_that("validating an in tree", {
  graph <- igraph::make_tree(3, children = 2, mode = "in")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must be an out tree")
})


test_that("validating graph without monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have a vertex attribute 'mono'")
})


test_that("validating graph with NA in monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have no NA in vertex attribute 'mono'")
})


test_that("validating graph without linkage attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have an edge attribute 'linkage'")
})
