test_that("low-level graph pipeline matches validated construction", {
  structures <- c(
    first = o_glycan_core_1(),
    second = n_glycan_core(),
    duplicate = o_glycan_core_1()
  )
  graphs <- get_structure_graphs(structures)

  graphs <- purrr::map(graphs, validate_glycan_graph)
  graphs <- purrr::map(graphs, canonicalize_glycan_graph)
  validate_glycan_graph_vector(graphs)
  iupacs <- purrr::map_chr(graphs, graph_to_iupac)

  unique <- !duplicated(unname(iupacs))
  unique_graphs <- graphs[unique]
  names(unique_graphs) <- unname(iupacs[unique])
  result <- new_glycan_structure(iupacs, unique_graphs)

  expected <- as_glycan_structure(unname(graphs))
  expect_equal(unname(structure_to_iupac(result)), structure_to_iupac(expected))
  expect_equal(names(result), names(structures))
  expect_length(attr(result, "graphs"), 2)
})


test_that("validate_glycan_graph returns valid input unchanged", {
  graph <- get_structure_graphs(o_glycan_core_1(), return_list = FALSE)

  expect_identical(validate_glycan_graph(graph), graph)
})


test_that("legacy scalar graph validation remains available", {
  graph <- get_structure_graphs(o_glycan_core_1(), return_list = FALSE)

  expect_identical(validate_single_glycan_structure(graph), graph)
})


test_that("validate_glycan_graph reports invalid input", {
  graph <- igraph::as_undirected(
    get_structure_graphs(o_glycan_core_1(), return_list = FALSE)
  )

  expect_snapshot(
    error = TRUE,
    validate_glycan_graph(graph)
  )
})


test_that("canonicalize_glycan_graph restores IUPAC ordering", {
  expected <- get_structure_graphs(n_glycan_core(), return_list = FALSE)
  scrambled <- igraph::permute(expected, rev(seq_len(igraph::vcount(expected))))

  result <- canonicalize_glycan_graph(scrambled)

  expect_equal(
    igraph::V(result)$name,
    as.character(seq_len(igraph::vcount(result)))
  )
  expect_equal(igraph::V(result)$mono, igraph::V(expected)$mono)
  expect_equal(igraph::E(result)$linkage, igraph::E(expected)$linkage)
})


test_that("validate_glycan_graph_vector rejects mixed monosaccharide types", {
  concrete <- get_structure_graphs(o_glycan_core_1(), return_list = FALSE)
  generic <- get_structure_graphs(
    convert_to_generic(o_glycan_core_1()),
    return_list = FALSE
  )

  expect_snapshot(
    error = TRUE,
    validate_glycan_graph_vector(list(concrete, generic))
  )
})


test_that("graph_to_iupac generates one string from one graph", {
  graph <- get_structure_graphs(n_glycan_core(), return_list = FALSE)

  expect_identical(
    graph_to_iupac(graph),
    unname(structure_to_iupac(n_glycan_core()))
  )
  expect_snapshot(
    error = TRUE,
    graph_to_iupac(list(graph))
  )
})


test_that("low-level graph pipeline supports annotated floating parts", {
  graph <- igraph::make_empty_graph(3, directed = TRUE)
  graph <- igraph::add_edges(graph, c(2, 3))
  igraph::V(graph)$mono <- c("Neu5Ac", "GalNAc", "Gal")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-3"
  graph$anomer <- "a1"
  graph$floating_parts <- list(
    list(root = 1L, linkage = "a2-6", parents = c(2L, 3L))
  )

  expect_identical(validate_glycan_graph(graph), graph)

  result <- canonicalize_glycan_graph(graph)

  expect_equal(igraph::V(result)$mono, c("Neu5Ac", "Gal", "GalNAc"))
  expect_equal(
    result$floating_parts,
    list(list(root = 1L, linkage = "a2-6", parents = c(2L, 3L)))
  )
  expect_identical(
    graph_to_iupac(result),
    "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
  )
})


test_that("unrestricted parts resolve against a one-node main tree", {
  graph <- igraph::make_empty_graph(2, directed = TRUE)
  igraph::V(graph)$mono <- c("Gal", "Neu5Ac")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- character()
  graph$anomer <- "a1"
  graph$floating_parts <- list(
    list(root = 2L, linkage = "a2-3", parents = integer())
  )

  result <- canonicalize_glycan_graph(validate_glycan_graph(graph))

  expect_identical(graph_to_iupac(result), "Neu5Ac(a2-3)Gal(a1-")
  expect_null(igraph::graph_attr(result, "floating_parts"))
})


test_that("unannotated forests remain invalid", {
  graph <- igraph::make_empty_graph(2, directed = TRUE)
  igraph::V(graph)$mono <- c("Gal", "Neu5Ac")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- character()
  graph$anomer <- "a1"

  expect_snapshot(
    error = TRUE,
    validate_glycan_graph(graph)
  )
})


test_that("new_glycan_structure checks graph lookup integrity", {
  graph <- get_structure_graphs(o_glycan_core_1(), return_list = FALSE)
  iupac <- graph_to_iupac(graph)

  result <- new_glycan_structure(
    c(first = iupac, duplicate = iupac),
    stats::setNames(list(graph), iupac)
  )

  expect_equal(names(result), c("first", "duplicate"))
  expect_equal(unname(structure_to_iupac(result)), rep(iupac, 2))

  partially_named <- rep(iupac, 2)
  attr(partially_named, "names") <- c("first", NA_character_)
  partially_named_result <- new_glycan_structure(
    partially_named,
    stats::setNames(list(graph), iupac)
  )
  expect_identical(names(partially_named_result), c("first", NA_character_))

  expect_snapshot(
    error = TRUE,
    new_glycan_structure(iupac, list(graph))
  )
  expect_snapshot(
    error = TRUE,
    new_glycan_structure("missing-key", stats::setNames(list(graph), iupac))
  )
})
