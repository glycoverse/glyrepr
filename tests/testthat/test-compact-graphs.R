test_that("native graph modes retain identity, attributes and serialization", {
  g <- .parse_iupac_condensed_single("Gal(b1-4)[Fuc(a1-3)]GlcNAc")
  g <- igraph::set_vertex_attr(g, "source", value = list("gal", 2L, c(3, 4)))
  g <- igraph::set_edge_attr(g, "weight", value = c(11, 22))
  g$origin <- list(format = "test")
  g <- igraph::permute(g, rev(seq_len(igraph::vcount(g))))
  expected <- testthat::with_mocked_bindings(
    canonicalize_graph_with_iupac(g),
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  result <- canonicalize_graph_with_iupac(g)
  expect_identical(validate_glycan_graph(g), g)
  expect_identical(result$iupac, expected$iupac)
  expect_identical(
    igraph::vertex_attr(result$graph),
    igraph::vertex_attr(expected$graph)
  )
  expect_identical(
    igraph::edge_attr(result$graph),
    igraph::edge_attr(expected$graph)
  )
  expect_identical(
    igraph::graph_attr(result$graph),
    igraph::graph_attr(expected$graph)
  )
  expect_identical(graph_to_iupac(result$graph), expected$iupac)
})

test_that("native transforms match the reference for floating and repeated inputs", {
  x <- as_glycan_structure(c(
    ordinary = "Gal(b1-4)[Fuc(a1-3)]GlcNAc-ol(?1-",
    floating = "{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc(?1-",
    substituent = "{3S|1,2}Gal(b1-4)Glc(?1-",
    modified = "Gal3S(b1-4)Glc(?1-",
    duplicate = "Gal(b1-4)[Fuc(a1-3)]GlcNAc-ol(?1-",
    missing = NA
  ))
  for (f in list(
    convert_to_generic,
    remove_linkages,
    remove_substituents,
    fill_anomer_pos
  )) {
    expected <- testthat::with_mocked_bindings(
      f(x),
      .compact_graph_results = function(graphs, ...) {
        vector("list", length(graphs))
      }
    )
    result <- f(x)
    expect_identical(as.character(result), as.character(expected))
    expect_identical(names(result), names(x))
    expect_identical(is.na(result), is.na(x))
    expect_equal(structure_nodes(result), structure_nodes(expected))
    expect_equal(structure_edges(result), structure_edges(expected))
    expect_equal(
      structure_floating_parts(result),
      structure_floating_parts(expected)
    )
    expect_equal(
      structure_floating_substituents(result),
      structure_floating_substituents(expected)
    )
  }
})

test_that("singleton canonicalization preserves attributes on original edges", {
  g <- igraph::make_graph(c(2, 3), n = 3, directed = TRUE)
  igraph::V(g)$mono <- c("Neu5Ac", "Glc", "Gal")
  igraph::V(g)$sub <- ""
  igraph::V(g)$source_id <- seq_len(3)
  igraph::E(g)$linkage <- "b1-4"
  igraph::E(g)$weight <- 42
  g$anomer <- "?1"
  g$floating_parts <- list(list(root = 1L, linkage = "a2-6", parents = 3L))
  expected <- testthat::with_mocked_bindings(
    canonicalize_glycan_graph(g),
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  result <- canonicalize_glycan_graph(g)
  expect_identical(igraph::as_edgelist(result), igraph::as_edgelist(expected))
  expect_identical(igraph::vertex_attr(result), igraph::vertex_attr(expected))
  expect_identical(igraph::edge_attr(result), igraph::edge_attr(expected))
  expect_identical(igraph::graph_attr(result), igraph::graph_attr(expected))
})

test_that("native serialization leaves singleton metadata unresolved", {
  g <- igraph::make_empty_graph(n = 2, directed = TRUE)
  igraph::V(g)$mono <- c("Gal", "Neu5Ac")
  igraph::V(g)$sub <- ""
  igraph::E(g)$linkage <- character()
  g$anomer <- "a1"
  g$floating_parts <- list(list(root = 2L, linkage = "a2-3", parents = 1L))
  expected <- testthat::with_mocked_bindings(
    graph_to_iupac(g),
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  expect_identical(graph_to_iupac(g), expected)
  expect_identical(validate_glycan_graph(g), g)
  expect_identical(
    graph_to_iupac(canonicalize_glycan_graph(g)),
    "Neu5Ac(a2-3)Gal(a1-"
  )
})

test_that("floating enumeration preserves candidate order and source IDs", {
  x <- as_glycan_structure("{3S|2,3}{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc(?1-")
  g <- as.list(x)[[1]]
  igraph::V(g)$source_id <- seq_len(igraph::vcount(g))
  igraph::E(g)$weight <- 12
  expected <- testthat::with_mocked_bindings(
    enumerate_floating_graph_localizations(g),
    .compact_enumerate_localizations = function(...) NULL,
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  result <- enumerate_floating_graph_localizations(g)
  expect_identical(result$assignments, expected$assignments)
  expect_identical(result$variant_id, expected$variant_id)
  for (i in seq_len(nrow(result))) {
    expect_identical(
      igraph::as_edgelist(result$graph[[i]]),
      igraph::as_edgelist(expected$graph[[i]])
    )
    expect_identical(
      igraph::vertex_attr(result$graph[[i]]),
      igraph::vertex_attr(expected$graph[[i]])
    )
    expect_identical(
      igraph::edge_attr(result$graph[[i]]),
      igraph::edge_attr(expected$graph[[i]])
    )
  }
})

test_that("valid native pipelines avoid the reference graph traversals", {
  x <- as_glycan_structure(c("Gal(b1-4)Glc(?1-", "{3S|1,2}Gal(b1-4)Glc(?1-"))
  g <- as.list(x)[[1]]
  testthat::local_mocked_bindings(
    build_seq_cache = function(...) stop("unexpected reference traversal"),
    validate_floating_graph_shape = function(...) {
      stop("unexpected reference validation")
    }
  )
  expect_identical(validate_glycan_graph(g), g)
  expect_identical(graph_to_iupac(g), unname(as.character(x)[1]))
  expect_identical(
    as.character(as_glycan_structure(as.list(x))),
    as.character(x)
  )
  expect_length(convert_to_generic(x), 2)
  expect_length(remove_linkages(x), 2)
  expect_length(remove_substituents(x), 2)
  expect_length(fill_anomer_pos(remove_linkages(x)), 2)
  expect_equal(nrow(enumerate_floating_localizations(x[2])), 2)
})

test_that("serialization preserves supplied candidate-parent order", {
  g <- as.list(as_glycan_structure("{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc(?1-"))[[1]]
  parts <- g$floating_parts
  parts[[1]]$parents <- rev(parts[[1]]$parents)
  g$floating_parts <- parts
  expected <- testthat::with_mocked_bindings(
    graph_to_iupac(g),
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  expect_identical(graph_to_iupac(g), expected)
  expect_identical(validate_glycan_graph(g), g)
})

test_that("table construction uses arrays and accepts omitted component nodes", {
  x <- as_glycan_structure(c(
    first = "{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc(?1-",
    missing = NA
  ))
  parts <- structure_floating_parts(x)
  parts$nodes <- NULL
  builder <- build_structure_graph_from_table_rows
  testthat::local_mocked_bindings(
    build_structure_graph_from_table_rows = function(..., .arrays = FALSE) {
      if (!.arrays) {
        stop("unexpected intermediate graph")
      }
      builder(..., .arrays = TRUE)
    }
  )
  result <- structure_from_tibbles(
    structure_nodes(x),
    structure_edges(x),
    get_anomer(x),
    floating_parts = parts,
    alditols = get_alditol(x)
  )
  expect_identical(as.character(result), as.character(x))
  expect_equal(structure_floating_parts(result), structure_floating_parts(x))
})

test_that("partial native localization preserves remaining metadata order", {
  x <- as_glycan_structure(
    "{Fuc(a1-2)|3,4}{Neu5Ac(a2-6)|3,4}Gal(b1-3)GalNAc(a1-"
  )
  g <- as.list(x)[[1]]
  parts <- g$floating_parts
  parts[[1]]$parents <- rev(parts[[1]]$parents)
  g$floating_parts <- parts
  assignments <- tibble::tibble(glycan_id = 1L, part_id = 2L, parent_node = 3L)
  expected <- testthat::with_mocked_bindings(
    localize_floating_parts(g, assignments),
    .compact_localize_parts = function(...) NULL,
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  result <- localize_floating_parts(g, assignments)
  expect_identical(igraph::as_edgelist(result), igraph::as_edgelist(expected))
  expect_identical(igraph::graph_attr(result), igraph::graph_attr(expected))
})

test_that("partial localization retains invalid graph diagnostics", {
  g <- as.list(as_glycan_structure("{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc(?1-"))[[1]]
  igraph::V(g)$mono[1] <- "unknown"
  assignments <- tibble::tibble(glycan_id = 1L, part_id = 1L, parent_node = 2L)
  condition <- function() {
    tryCatch(localize_floating_parts(g, assignments), error = identity)
  }
  expected <- testthat::with_mocked_bindings(
    condition(),
    .compact_localize_parts = function(...) NULL,
    .compact_graph_results = function(graphs, ...) {
      vector("list", length(graphs))
    }
  )
  result <- condition()
  expect_s3_class(result, "error")
  expect_identical(conditionMessage(result), conditionMessage(expected))
})
