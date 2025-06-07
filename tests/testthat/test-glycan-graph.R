test_that("glycan graph class", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")

  glycan <- new_glycan_graph(graph)

  expect_s3_class(glycan, c("glycan_graph", "igraph"))
})


test_that("validating undirected graphs", {
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must be directed")
})


test_that("validating an in tree", {
  graph <- igraph::make_tree(3, children = 2, mode = "in")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must be an out tree")
})


test_that("validating graph without monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have a vertex attribute 'mono'")
})


test_that("validating graph without substituent attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have a vertex attribute 'sub'")
})


test_that("validating graph with NA in monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have no NA in vertex attribute 'mono'")
})


test_that("validating graph with NA in substitude attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- c("", NA, "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have no NA in vertex attribute 'sub'")
})


patrick::with_parameters_test_that("valid substituents", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- sub
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "b1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_no_error(validate_glycan_graph(glycan))
}, sub = c("6S", "9Ac", "2P", "?S"))


test_that("validating graph without linkage attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have an edge attribute 'linkage'")
})


test_that("validating one non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Fuc", "Bad")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan))
  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_equal(err$message, "Unknown monosaccharide: Bad")
  expect_equal(err$monos, "Bad")
})


test_that("validating two non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Bad1", "Bad2")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_equal(err$message, "Unknown monosaccharide: Bad1, Bad2")
  expect_equal(err$monos, c("Bad1", "Bad2"))
})


test_that("validating duplicated non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Bad", "Bad")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_equal(err$message, "Unknown monosaccharide: Bad")
  expect_equal(err$monos, "Bad")
})


test_that("validating bad subtituent", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- c("", "6S", "Bad")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan))
  err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  expect_equal(err$message, "Unknown substituent: Bad")
  expect_equal(err$subs, "Bad")
})


patrick::with_parameters_test_that("validating bad linkage", {
    graph <- igraph::make_graph(~ 1-+2, 2-+3)
    igraph::V(graph)$mono <- c("Hex", "Fuc", "Hex")
    igraph::V(graph)$sub <- ""
    igraph::E(graph)$linkage <- bad_linkage
    graph$anomer <- "a1"
    graph$alditol <- FALSE

    glycan <- new_glycan_graph(graph)

    expect_error(validate_glycan_graph(glycan))
    err <- rlang::catch_cnd(validate_glycan_graph(glycan))
  },
  bad_linkage = c("1-4", "c1-4", "b1", "abc", ""),
  .test_name = bad_linkage
)


test_that("validating NA linkages", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", NA)
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan))
})


test_that("validating mixed generic and concrete monosaccharides", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "GlcNAc", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Monosaccharides must be either all generic or all concrete")
})


test_that("missing anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have a graph attribute 'anomer'")
})


test_that("invalid anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a"
  graph$alditol <- FALSE

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Invalid anomer: a")
})


test_that("missing alditol attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph must have a graph attribute 'alditol'")
})


test_that("invalid alditol attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- "True"

  glycan <- new_glycan_graph(graph)

  expect_error(validate_glycan_graph(glycan), "Glycan graph attribute 'alditol' must be logical")
})
