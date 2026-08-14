.make_symmetric_floating_graph <- function(
  edge_order = c(1, 2, 1, 3),
  parents = c(1L, 2L)
) {
  graph <- igraph::make_empty_graph(4, directed = TRUE)
  graph <- igraph::add_edges(graph, edge_order)
  igraph::V(graph)$name <- as.character(seq_len(4))
  igraph::V(graph)$mono <- c("Glc", "Gal", "Gal", "Neu5Ac")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("a1-?", "a1-?")
  graph$anomer <- "a1"
  graph$floating_parts <- list(
    list(root = 4L, linkage = "a2-3", parents = parents)
  )
  graph
}


.floating_symmetry_parent_positions <- function(graph) {
  info <- floating_graph_info(graph)
  order <- floating_symmetry_main_order(graph, info)
  canonical_names <- igraph::V(graph)$name[order$vertices]
  canonical_index <- stats::setNames(
    seq_along(canonical_names),
    canonical_names
  )

  purrr::map(
    info$parts,
    function(part) {
      parent_names <- igraph::V(graph)$name[part$parents]
      sort(as.integer(unname(canonical_index[parent_names])))
    }
  )
}


test_that("symmetric edge permutations normalize explicit parents", {
  graph_a <- .make_symmetric_floating_graph(c(1, 2, 1, 3))
  graph_b <- .make_symmetric_floating_graph(c(1, 3, 1, 2))

  expect_equal(
    .floating_symmetry_parent_positions(graph_a),
    list(c(1L, 3L))
  )
  expect_equal(
    .floating_symmetry_parent_positions(graph_b),
    list(c(1L, 3L))
  )

  glycan_a <- glycan_structure(graph_a)
  glycan_b <- glycan_structure(graph_b)
  expected <- "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-"

  expect_identical(structure_to_iupac(glycan_a), expected)
  expect_identical(structure_to_iupac(glycan_b), expected)
  expect_true(unname(glycan_a == glycan_b))
})


test_that("symmetric vertex permutations normalize explicit parents", {
  graph_a <- .make_symmetric_floating_graph()

  graph_b <- igraph::make_empty_graph(4, directed = TRUE)
  graph_b <- igraph::add_edges(graph_b, c(3, 4, 3, 1))
  igraph::V(graph_b)$name <- c("gal-a", "floating", "root", "gal-b")
  igraph::V(graph_b)$mono <- c("Gal", "Neu5Ac", "Glc", "Gal")
  igraph::V(graph_b)$sub <- ""
  igraph::E(graph_b)$linkage <- c("a1-?", "a1-?")
  graph_b$anomer <- "a1"
  graph_b$floating_parts <- list(
    list(root = 2L, linkage = "a2-3", parents = c(1L, 3L))
  )

  expect_equal(
    .floating_symmetry_parent_positions(graph_a),
    .floating_symmetry_parent_positions(graph_b)
  )

  glycan_a <- glycan_structure(graph_a)
  glycan_b <- glycan_structure(graph_b)

  expect_identical(structure_to_iupac(glycan_a), structure_to_iupac(glycan_b))
  expect_true(unname(glycan_a == glycan_b))
})


test_that("floating canonicalization does not depend on input vertex names", {
  graph <- .make_symmetric_floating_graph()
  igraph::V(graph)$name <- rep("duplicate", igraph::vcount(graph))

  glycan <- glycan_structure(graph)

  expect_identical(
    structure_to_iupac(glycan),
    "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-"
  )
})


test_that("nested symmetric branches use descendant parent constraints", {
  make_graph <- function(edge_order) {
    graph <- igraph::make_empty_graph(6, directed = TRUE)
    graph <- igraph::add_edges(graph, edge_order)
    igraph::V(graph)$name <- as.character(seq_len(6))
    igraph::V(graph)$mono <- c(
      "Glc",
      "GlcNAc",
      "GlcNAc",
      "Gal",
      "Gal",
      "Neu5Ac"
    )
    igraph::V(graph)$sub <- ""
    igraph::E(graph)$linkage <- c("b1-?", "b1-4", "b1-?", "b1-4")
    graph$anomer <- "a1"
    graph$floating_parts <- list(
      list(root = 6L, linkage = "a2-3", parents = c(1L, 4L))
    )
    graph
  }

  graph_a <- make_graph(c(1, 2, 2, 4, 1, 3, 3, 5))
  graph_b <- make_graph(c(1, 3, 3, 5, 1, 2, 2, 4))

  expect_equal(
    .floating_symmetry_parent_positions(graph_a),
    list(c(1L, 5L))
  )
  expect_equal(
    .floating_symmetry_parent_positions(graph_b),
    list(c(1L, 5L))
  )

  glycan_a <- glycan_structure(graph_a)
  glycan_b <- glycan_structure(graph_b)

  expect_identical(structure_to_iupac(glycan_a), structure_to_iupac(glycan_b))
  expect_true(unname(glycan_a == glycan_b))
})


test_that("multiple parent sets are canonicalized jointly", {
  make_graph <- function(same_arm, edge_order) {
    graph <- igraph::make_empty_graph(5, directed = TRUE)
    graph <- igraph::add_edges(graph, edge_order)
    igraph::V(graph)$name <- as.character(seq_len(5))
    igraph::V(graph)$mono <- c("Glc", "Gal", "Gal", "Fuc", "Neu5Ac")
    igraph::V(graph)$sub <- ""
    igraph::E(graph)$linkage <- c("a1-?", "a1-?")
    graph$anomer <- "a1"
    graph$floating_parts <- list(
      list(root = 4L, linkage = "a1-6", parents = c(1L, 2L)),
      list(
        root = 5L,
        linkage = "a2-3",
        parents = if (same_arm) c(1L, 2L) else c(1L, 3L)
      )
    )
    graph
  }

  same_a <- make_graph(TRUE, c(1, 2, 1, 3))
  same_b <- make_graph(TRUE, c(1, 3, 1, 2))
  split_a <- make_graph(FALSE, c(1, 2, 1, 3))
  split_b <- make_graph(FALSE, c(1, 3, 1, 2))

  expect_equal(
    .floating_symmetry_parent_positions(same_a),
    list(c(1L, 3L), c(1L, 3L))
  )
  expect_equal(
    .floating_symmetry_parent_positions(same_b),
    list(c(1L, 3L), c(1L, 3L))
  )
  expect_equal(
    .floating_symmetry_parent_positions(split_a),
    .floating_symmetry_parent_positions(split_b)
  )
  expect_false(
    identical(
      .floating_symmetry_parent_positions(same_a),
      .floating_symmetry_parent_positions(split_a)
    )
  )

  same_glycan_a <- glycan_structure(same_a)
  same_glycan_b <- glycan_structure(same_b)
  split_glycan_a <- glycan_structure(split_a)
  split_glycan_b <- glycan_structure(split_b)

  expect_true(unname(same_glycan_a == same_glycan_b))
  expect_true(unname(split_glycan_a == split_glycan_b))
  expect_false(unname(same_glycan_a == split_glycan_a))
})


test_that("cross-component domains ignore component and metadata order", {
  graph_a <- igraph::make_empty_graph(3, directed = TRUE)
  igraph::V(graph_a)$name <- c("fuc", "man", "glc")
  igraph::V(graph_a)$mono <- c("Fuc", "Man", "Glc")
  igraph::V(graph_a)$sub <- ""
  igraph::E(graph_a)$linkage <- character()
  graph_a$anomer <- "a1"
  graph_a$floating_parts <- list(
    list(root = 1L, linkage = "a1-2", parents = c(2L, 3L)),
    list(root = 2L, linkage = "a1-3", parents = c(1L, 3L))
  )
  graph_a$floating_substituents <- list(
    list(substituent = "6S", parents = c(1L, 3L))
  )

  graph_b <- igraph::make_empty_graph(3, directed = TRUE)
  igraph::V(graph_b)$name <- c("glc", "man", "fuc")
  igraph::V(graph_b)$mono <- c("Glc", "Man", "Fuc")
  igraph::V(graph_b)$sub <- ""
  igraph::E(graph_b)$linkage <- character()
  graph_b$anomer <- "a1"
  graph_b$floating_parts <- list(
    list(root = 2L, linkage = "a1-3", parents = c(3L, 1L)),
    list(root = 3L, linkage = "a1-2", parents = c(2L, 1L))
  )
  graph_b$floating_substituents <- list(
    list(substituent = "6S", parents = c(3L, 1L))
  )

  glycan_a <- glycan_structure(graph_a)
  glycan_b <- glycan_structure(graph_b)

  expect_identical(structure_to_iupac(glycan_a), structure_to_iupac(glycan_b))
  expect_true(unname(glycan_a == glycan_b))
})
