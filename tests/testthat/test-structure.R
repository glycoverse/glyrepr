# Tests for glycan structure functions

good_glycan_structure <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  graph
}


# Tests for glycan_structure --------------------------------------------------

test_that("glycan_structure works", {
  glycan <- glycan_structure(good_glycan_structure())
  expect_s3_class(glycan, c("glycan_structure", "igraph"))
})


test_that("glycan_structure fails for invalid graphs", {
  bad_graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+1)
  expect_error(glycan_structure(bad_graph))
})


test_that("vertex names are added if missing", {
  graph <- good_glycan_structure()
  graph <- igraph::delete_vertex_attr(graph, "name")

  glycan <- glycan_structure(graph)

  expect_true("name" %in% igraph::vertex_attr_names(glycan))
})


# Tests for new_glycan_structure -------------------------------------------------

test_that("new_glycan_structure creates correct class", {
  graph <- igraph::make_empty_graph()
  glycan <- new_glycan_structure(graph)
  expect_s3_class(glycan, c("glycan_structure", "igraph"))
})


test_that("new_glycan_structure requires igraph input", {
  expect_error(new_glycan_structure("not a graph"))
})


# Tests for validate_glycan_structure --------------------------------------------

test_that("validate_glycan_structure accepts valid graphs", {
  graph <- good_glycan_structure()
  glycan <- new_glycan_structure(graph)
  expect_no_error(validate_glycan_structure(glycan))
})


test_that("validate_glycan_structure rejects invalid graphs", {
  # Test undirected graph
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("Glc", "Gal")
  igraph::V(graph)$sub <- c("", "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_glycan_structure(graph)
  expect_error(validate_glycan_structure(glycan), "directed")
})


# Tests for ensure_name_vertex_attr ------------------------------------------

test_that("ensure_name_vertex_attr adds names when missing", {
  graph <- igraph::make_graph(~ 1-+2)
  graph <- igraph::delete_vertex_attr(graph, "name")
  result <- ensure_name_vertex_attr(graph)
  expect_true("name" %in% igraph::vertex_attr_names(result))
})


test_that("ensure_name_vertex_attr preserves existing names", {
  graph <- igraph::make_graph(~ A-+B)
  result <- ensure_name_vertex_attr(graph)
  expect_equal(igraph::V(result)$name, c("A", "B"))
})


test_that("glycan structure class", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")

  glycan <- new_glycan_structure(graph)

  expect_s3_class(glycan, c("glycan_structure", "igraph"))
})


test_that("validating undirected graphs", {
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must be directed")
})


test_that("validating an in tree", {
  graph <- igraph::make_tree(3, children = 2, mode = "in")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must be an out tree")
})


test_that("validating graph without monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have a vertex attribute 'mono'")
})


test_that("validating graph without substituent attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have a vertex attribute 'sub'")
})


test_that("validating graph with NA in monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have no NA in vertex attribute 'mono'")
})


test_that("validating graph with NA in substitude attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- c("", NA, "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have no NA in vertex attribute 'sub'")
})


patrick::with_parameters_test_that("valid substituents", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- sub
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "b1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_no_error(validate_glycan_structure(glycan))
}, sub = c("6S", "9Ac", "2P", "?S"))


test_that("validating graph without linkage attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have an edge attribute 'linkage'")
})


test_that("validating one non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Fuc", "Bad")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan))
  err <- rlang::catch_cnd(validate_glycan_structure(glycan))
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

  glycan <- new_glycan_structure(graph)

  err <- rlang::catch_cnd(validate_glycan_structure(glycan))
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

  glycan <- new_glycan_structure(graph)

  err <- rlang::catch_cnd(validate_glycan_structure(glycan))
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

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan))
  err <- rlang::catch_cnd(validate_glycan_structure(glycan))
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

    glycan <- new_glycan_structure(graph)

    expect_error(validate_glycan_structure(glycan))
    err <- rlang::catch_cnd(validate_glycan_structure(glycan))
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

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan))
})


test_that("validating mixed generic and concrete monosaccharides", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "GlcNAc", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Monosaccharides must be either all generic or all concrete")
})


test_that("missing anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have a graph attribute 'anomer'")
})


test_that("invalid anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a"
  graph$alditol <- FALSE

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Invalid anomer: a")
})


test_that("missing alditol attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure must have a graph attribute 'alditol'")
})


test_that("invalid alditol attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- "True"

  glycan <- new_glycan_structure(graph)

  expect_error(validate_glycan_structure(glycan), "Glycan structure attribute 'alditol' must be logical")
})
