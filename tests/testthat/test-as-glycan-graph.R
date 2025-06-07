good_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  graph
}


test_that("as_glycan_graph works", {
  glycan <- as_glycan_graph(good_glycan_graph())
  expect_s3_class(glycan, c("glycan_graph", "igraph"))
})


test_that("as_glycan_graph fails for invalid graphs", {
  bad_graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+1)
  expect_error(as_glycan_graph(bad_graph))
})


test_that("vertex names are added if missing", {
  graph <- good_glycan_graph()
  graph <- igraph::delete_vertex_attr(graph, "name")

  glycan <- as_glycan_graph(graph)

  expect_true("name" %in% igraph::vertex_attr_names(glycan))
})
