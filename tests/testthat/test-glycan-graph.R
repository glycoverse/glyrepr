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


test_that("validating one non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Fuc", "Bad")
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan))
  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_s3_class(err, "error_bad_mono")
  expect_equal(err$message, "Unknown monosaccharide: Bad")
  expect_equal(err$monos, "Bad")
})


test_that("validating two non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Bad1", "Bad2")
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_equal(err$message, "Unknown monosaccharide: Bad1, Bad2")
  expect_equal(err$monos, c("Bad1", "Bad2"))
})


test_that("validating duplicated non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Bad", "Bad")
  igraph::E(graph)$linkage <- "b1-4"

  glycan <- new_glycan_graph(graph)

  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_equal(err$message, "Unknown monosaccharide: Bad")
  expect_equal(err$monos, "Bad")
})


patrick::with_parameters_test_that("validating bad linkage", {
    graph <- igraph::make_graph(~ 1-+2, 2-+3)
    igraph::V(graph)$mono <- c("Hex", "Fuc", "Hex")
    igraph::E(graph)$linkage <- bad_linkage

    glycan <- new_glycan_graph(graph)

    expect_error(validate_glycan_graph(glycan))
    err <- rlang::catch_cnd(validate_glycan_graph(glycan))
    expect_s3_class(err, "error_bad_linkage")
  },
  bad_linkage = c("1-4", "c1-4", "b1", "abc", ""),
  .test_name = bad_linkage
)
