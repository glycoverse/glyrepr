good_dn_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Glc")
  igraph::V(graph)$sub <- c("", NA, "")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA)
  graph$anomer <- "a1"
  graph
}


good_ne_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph
}


test_that("as_dn_glycan_graph works", {
  glycan <- as_dn_glycan_graph(good_dn_graph())
  expect_s3_class(glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
})


test_that("as_ne_glycan_graph works", {
  glycan <- as_ne_glycan_graph(good_ne_graph())
  expect_s3_class(glycan, c("ne_glycan_graph", "glycan_graph", "igraph"))
})


test_that("as_glycan_graph works for DN graphs", {
  glycan <- as_glycan_graph(good_dn_graph(), type = "dn")
  expect_s3_class(glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
})


test_that("as_glycan_graph works for NE graphs", {
  glycan <- as_glycan_graph(good_ne_graph(), type = "ne")
  expect_s3_class(glycan, c("ne_glycan_graph", "glycan_graph", "igraph"))
})


test_that("as_glycan_graph works for DN graphs (auto)", {
  dn_glycan <- as_glycan_graph(good_dn_graph())
  expect_s3_class(dn_glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
})


test_that("as_glycan_graph works for NE graphs (auto)", {
  ne_glycan <- as_glycan_graph(good_ne_graph())
  expect_s3_class(ne_glycan, c("ne_glycan_graph", "glycan_graph", "igraph"))
})


test_that("as_glycan_graph fails for invalid graphs", {
  bad_graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+1)
  expect_error(as_glycan_graph(bad_graph), "Could not infer glycan graph type")
})


test_that("vertex names are added if missing", {
  graph <- good_ne_graph()
  graph <- igraph::delete_vertex_attr(graph, "name")

  glycan <- as_ne_glycan_graph(graph)

  expect_true("name" %in% igraph::vertex_attr_names(glycan))
})
