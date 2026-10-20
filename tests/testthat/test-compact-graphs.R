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
      .compact_transform_graphs_native = function(...) NULL,
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

test_that("batch transforms preserve custom metadata and first collision graphs", {
  graphs <- as.list(as_glycan_structure(c(
    "Gal3S(b1-4)[Fuc(a1-3)]GlcNAc",
    "{3S|1,2}Gal(b1-4)Glc",
    "{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc",
    "Gal(b1-4)GlcNAc",
    "Glc(b1-4)GlcNAc"
  )))
  graphs[[1]] <- igraph::set_vertex_attr(
    graphs[[1]],
    "source",
    value = seq_len(igraph::vcount(graphs[[1]]))
  )
  graphs[[1]] <- igraph::set_edge_attr(graphs[[1]], "weight", value = c(4, 8))
  graphs[[1]] <- igraph::set_graph_attr(
    graphs[[1]],
    "origin",
    list(source = "test")
  )
  graphs[[4]] <- igraph::set_graph_attr(graphs[[4]], "origin", "first")
  graphs[[5]] <- igraph::set_graph_attr(graphs[[5]], "origin", "second")
  x <- as_glycan_structure(graphs)[c(1, 2, 3, 4, 5, 1, NA)]
  names(x) <- letters[seq_along(x)]
  describe <- function(g) {
    list(
      edges = igraph::as_edgelist(g, names = FALSE),
      vertex_attrs = igraph::vertex_attr(g),
      edge_attrs = igraph::edge_attr(g),
      graph_attrs = igraph::graph_attr(g)
    )
  }
  original <- lapply(attr(x, "graphs"), describe)
  for (f in list(remove_linkages, convert_to_generic, remove_substituents)) {
    expected <- with_mocked_bindings(
      f(x),
      .compact_transform_graphs_native = function(...) NULL
    )
    actual <- f(x)
    expect_identical(as.character(actual), as.character(expected))
    expect_identical(names(actual), names(x))
    expect_identical(
      lapply(attr(actual, "graphs"), describe),
      lapply(attr(expected, "graphs"), describe)
    )
    expect_identical(lapply(attr(x, "graphs"), describe), original)
    expect_identical(f(actual), actual)
  }
})

test_that("batch transforms retain graph behavior after serialization and mutation", {
  x <- as_glycan_structure(c(
    "Gal3S(b1-4)[Fuc(a1-3)]GlcNAc-ol",
    "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-",
    "{6S|1,2}Gal(b1-4)Glc",
    "Gal",
    NA
  ))
  # Materialize igraph's native indexes before transforming the same graph data.
  invisible(lapply(attr(x, "graphs"), igraph::degree))
  describe <- function(g) {
    list(
      edges = igraph::as_edgelist(g, names = FALSE),
      vertex_attrs = igraph::vertex_attr(g),
      edge_attrs = igraph::edge_attr(g),
      graph_attrs = igraph::graph_attr(g),
      degree = igraph::degree(g)
    )
  }
  for (f in list(remove_linkages, convert_to_generic, remove_substituents)) {
    expected <- with_mocked_bindings(
      f(x),
      .compact_transform_graphs_native = function(...) NULL
    )
    actual <- unserialize(serialize(f(x), NULL))
    expect_identical(as.character(actual), as.character(expected))
    expect_identical(
      lapply(attr(actual, "graphs"), describe),
      lapply(attr(expected, "graphs"), describe)
    )
    g <- attr(actual, "graphs")[[1]]
    changed <- igraph::set_vertex_attr(g, "mono", 1, "Man")
    changed <- igraph::add_vertices(changed, 1)
    expect_equal(igraph::vcount(changed), igraph::vcount(g) + 1)
    expect_identical(describe(g), describe(attr(expected, "graphs")[[1]]))
    expect_identical(f(actual), actual)
  }
})

test_that("batch transforms fall back for other graph schemas and attributed fields", {
  x <- as_glycan_structure("Gal3S(b1-4)Glc")
  for (f in list(remove_linkages, convert_to_generic, remove_substituents)) {
    expected <- f(x)
    actual <- with_mocked_bindings(f(x), .compact_graph_template = function() {
      NULL
    })
    expect_identical(as.character(actual), as.character(expected))
    expect_identical(
      structure_nodes(actual),
      structure_nodes(expected)
    )
    expect_identical(structure_edges(actual), structure_edges(expected))
    expect_identical(f(x[FALSE]), x[FALSE])
    expect_identical(f(x[NA_integer_]), x[NA_integer_])
  }
  g <- as.list(x)[[1]]
  igraph::V(g)$sub <- structure(c("3S", ""), source = "test")
  x <- new_glycan_structure(as.character(x), setNames(list(g), as.character(x)))
  for (f in list(remove_linkages, convert_to_generic, remove_substituents)) {
    expected <- with_mocked_bindings(
      f(x),
      .compact_transform_graphs_native = function(...) NULL
    )
    actual <- f(x)
    expect_identical(as.character(actual), as.character(expected))
    expect_identical(
      igraph::vertex_attr(as.list(actual)[[1]]),
      igraph::vertex_attr(as.list(expected)[[1]])
    )
  }
})

test_that("batch graph import rejects incompatible fields before native traversal", {
  x <- as_glycan_structure("{Neu5Ac(a2-6)|2,3}Gal(b1-4)Glc")
  original <- unclass(as.list(x)[[1]])
  bad_endpoint <- original
  bad_endpoint[[3]][[1]] <- 1000
  bad_root <- original
  bad_root[[9]][[2]]$floating_parts[[1]]$root <- NA_integer_
  bad_nodes <- original
  bad_nodes[[9]][[2]]$floating_parts[[1]]$nodes <- 1000L
  bad_linkage <- original
  bad_linkage[[9]][[4]]$linkage[[1]] <- NA_character_
  bad_flags <- original
  bad_flags[[9]][[1]] <- NULL
  graphs <- lapply(
    list(bad_endpoint, bad_root, bad_nodes, bad_linkage, bad_flags),
    structure,
    class = "igraph"
  )
  graphs <- c(
    graphs,
    list(structure(original, class = "igraph", source = "test"))
  )
  result <- .compact_transform_graphs_native(
    graphs,
    rep(as.character(x), length(graphs)),
    "remove_linkages",
    character(),
    character(),
    base::order,
    .compact_bliss_labels,
    FALSE,
    .compact_graph_template()
  )
  skip_if(is.null(result), "Installed igraph uses an unsupported graph format")
  expect_identical(result$supported, rep(FALSE, length(graphs)))
  expect_identical(result$graphs, rep(list(NULL), length(graphs)))
})

test_that("batch graph metadata remains independent when a result is edited", {
  x <- as_glycan_structure(c("Gal(b1-4)Glc", "Man(a1-3)Glc"))
  y <- convert_to_generic(x)
  graphs <- as.list(y)
  before <- lapply(graphs, igraph::graph_attr)
  changed <- igraph::set_graph_attr(graphs[[1]], "anomer", "b1")
  changed <- igraph::set_graph_attr(changed, "alditol", TRUE)
  expect_identical(igraph::graph_attr(changed, "anomer"), "b1")
  expect_identical(igraph::graph_attr(changed, "alditol"), TRUE)
  expect_identical(lapply(graphs, igraph::graph_attr), before)
  expect_identical(lapply(as.list(y), igraph::graph_attr), before)
})

test_that("batch transforms normalize singleton metadata on low-level inputs", {
  g <- as.list(as_glycan_structure("Gal(b1-4)Glc"))[[1]]
  g$floating_substituents <- list(list(substituent = "3S", parents = 1L))
  key <- graph_to_iupac(g)
  x <- new_glycan_structure(key, setNames(list(g), key))
  for (f in list(convert_to_generic, remove_linkages, remove_substituents)) {
    expected <- with_mocked_bindings(
      f(x),
      .compact_transform_graphs_native = function(...) NULL
    )
    actual <- f(x)
    expect_identical(as.character(actual), as.character(expected))
    expect_identical(structure_nodes(actual), structure_nodes(expected))
    expect_identical(structure_edges(actual), structure_edges(expected))
    expect_identical(
      structure_floating_substituents(actual),
      structure_floating_substituents(expected)
    )
  }
})

test_that("batch transforms retain locale-specific ordering between calls", {
  withr::local_locale(c(LC_COLLATE = "C"))
  for (locale in c("C", "C.UTF-8", "en_US.UTF-8", "C")) {
    if (!nzchar(suppressWarnings(Sys.setlocale("LC_COLLATE", locale)))) {
      next
    }
    x <- as_glycan_structure(c(
      "Galf(b1-?)[GlcNAc(b1-?)]Man",
      "Galf(b1-?)[GlcNAc(b1-?)]Gal",
      "{Neu5Ac(a2-?)|2,3}Gal(b1-4)Glc",
      "Gal3S(b1-?)[GlcNAc(b1-?)]Man"
    ))
    for (f in list(convert_to_generic, remove_linkages, remove_substituents)) {
      expected <- with_mocked_bindings(
        f(x),
        .compact_transform_graphs_native = function(...) NULL
      )
      actual <- f(x)
      expect_identical(as.character(actual), as.character(expected))
      expect_identical(structure_nodes(actual), structure_nodes(expected))
      expect_identical(structure_edges(actual), structure_edges(expected))
      expect_identical(
        structure_floating_parts(actual),
        structure_floating_parts(expected)
      )
    }
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
