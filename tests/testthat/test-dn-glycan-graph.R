test_that("dual-node glycan graph class", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")

  glycan <- new_ne_glycan_graph(graph)

  expect_s3_class(glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
})


test_that("check if it is a directed graph", {
  graph <- igraph::make_graph(~ 1--2, 2--3)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Glycan graph must be directed.")
})


test_that("check if it is an out tree", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+1)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Glycan graph must be an out tree.")
})


test_that("check missing 'type' vertex attribute", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Glycan graph must have vertex attributes 'type', 'mono' and 'linkage'.")
})


test_that("check missing 'mono' vertex attribute", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Glycan graph must have vertex attributes 'type', 'mono' and 'linkage'.")
})


test_that("check missing 'linkage' vertex attribute", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Glycan graph must have vertex attributes 'type', 'mono' and 'linkage'.")
})


test_that("check NA in 'type' attribute", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", NA, "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Glycan graph must have no NA in 'type' attribute.")
})


test_that("check bad 'type' attribute", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "bad")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "The 'type' of a node could either be 'mono' or 'linkage'.")
})


test_that("check if 'mono' and 'linkage' nodes are alternating 1", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "mono", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")
  igraph::V(graph)$linkage <- c(NA, NA, NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "The 'mono' and 'linkage' nodes must be alternating.")
})


test_that("check if 'mono' and 'linkage' nodes are alternating 2", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("linkage", "linkage", "linkage")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")
  igraph::V(graph)$linkage <- c(NA, NA, NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "The 'mono' and 'linkage' nodes must be alternating.")
})


test_that("check if 'mono' and 'linkage' nodes are alternating 3", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "linkage")
  igraph::V(graph)$mono <- c("Glc", NA, NA)
  igraph::V(graph)$linkage <- c(NA, NA, NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "The 'mono' and 'linkage' nodes must be alternating.")
})


test_that("check NA in 'mono' attribute for 'mono' nodes", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, NA)
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Mono nodes must have no NA in 'mono' attribute.")
})


test_that("check unknown monosaccharides", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Xxx")
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Unknown monosaccharide: Xxx")
})


test_that("check mixed use of generic and concrete monosaccharides", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Hex")
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Mono nodes must not mix generic and concrete monosaccharides.")
})


test_that("check invalid linkages", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal")
  igraph::V(graph)$linkage <- c(NA, "a1", NA)

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan), "Invalid linkage: a1")
})


test_that("check 0 out degree linkage nodes", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4)
  igraph::V(graph)$type <- c("mono", "linkage", "mono", "linkage")
  igraph::V(graph)$mono <- c("Glc", NA, "Gal", NA)
  igraph::V(graph)$linkage <- c(NA, "a1-2", NA, "a1-2")

  glycan <- new_dn_glycan_graph(graph)

  expect_error(validate_dn_glycan_graph(glycan))
})
