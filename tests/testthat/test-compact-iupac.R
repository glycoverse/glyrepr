test_that("character recovery preserves repeated failures and canonical duplicates", {
  x <- c(
    modified = "Gal6S",
    bad = "not-a-structure",
    missing = NA_character_,
    bad_again = "not-a-structure",
    floating = "{6S|1}Gal",
    alditol = "Gal-ol",
    unknown = "Hex(??-"
  )
  expect_snapshot(result <- as_glycan_structure(x, on_failure = "na"))
  expect_identical(
    as.character(result),
    c(
      modified = "Gal6S(?1-",
      bad = NA_character_,
      missing = NA_character_,
      bad_again = NA_character_,
      floating = "Gal6S(?1-",
      alditol = "Gal-ol(?1-",
      unknown = "Hex(??-"
    )
  )
  expect_identical(length(attr(result, "graphs")), 3L)
  expect_identical(names(result), names(x))
})

test_that("strict character errors retain parse-before-validation precedence", {
  expect_snapshot(
    error = TRUE,
    as_glycan_structure(c(
      "Gal(b1-4)[Man(a1-4)]Glc",
      "not-a-structure"
    ))
  )
})

test_that("native results preserve floating metadata through related APIs", {
  x <- c(
    "Gal6S(b1-4)GlcNAc-ol",
    "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-",
    "{6S|1,2}Gal(b1-4)Glc",
    NA_character_
  )
  reference <- suppressWarnings(as_glycan_structure(
    lapply(x, function(s) {
      if (is.na(s)) {
        return(NA)
      }
      .parse_iupac_condensed_single(s)
    }),
    on_failure = "na"
  ))
  for (policy in c("error", "na")) {
    result <- as_glycan_structure(x, on_failure = policy)
    expect_identical(as.character(result), as.character(reference))
    for (accessor in list(
      structure_nodes,
      structure_edges,
      structure_floating_parts,
      structure_floating_substituents,
      structure_component_membership,
      structure_floating_candidates,
      get_alditol,
      get_anomer,
      get_mono_type,
      get_structure_level,
      has_linkages,
      count_mono
    )) {
      expect_equal(accessor(result), accessor(reference))
    }
    expect_identical(
      as.character(convert_to_generic(result)),
      as.character(convert_to_generic(reference))
    )
    expect_identical(
      as.character(remove_linkages(result)),
      as.character(remove_linkages(reference))
    )
  }
})

test_that("graph conversion retains arbitrary graph vertex and edge attributes", {
  graph <- .parse_iupac_condensed_single("Gal(b1-4)Glc")
  graph <- igraph::set_graph_attr(graph, "source", list(id = "example"))
  graph <- igraph::set_vertex_attr(
    graph,
    "label",
    value = c("reducing", "terminal")
  )
  graph <- igraph::set_edge_attr(graph, "weight", value = 2.5)
  for (policy in c("error", "na")) {
    result <- get_structure_graphs(as_glycan_structure(
      graph,
      on_failure = policy
    ))
    expect_identical(result$source, list(id = "example"))
    expect_identical(igraph::V(result)$label, c("terminal", "reducing"))
    expect_identical(igraph::E(result)$weight, 2.5)
  }
})

test_that("linkage scanning preserves unknown alternatives and floating linkages", {
  x <- c(
    "Gal(b1-3/?/6)Glc",
    "Gal(?-?)Glc",
    "{Gal(b1-3/?/6)}Glc",
    "{Gal(?-?)}Glc"
  )
  expect_identical(
    as.character(as_glycan_structure(x)),
    c(
      "Gal(b1-?)Glc(?1-",
      "Gal(??-?)Glc(?1-",
      "Gal(b1-?)Glc(?1-",
      "Gal(??-?)Glc(?1-"
    )
  )
  bad <- c("Gal(b1-?/?/6)Glc", "Gal(b1-3/)Glc", "Gal(b1-03)Glc")
  expect_snapshot(result <- as_glycan_structure(bad, on_failure = "na"))
  expect_identical(is.na(result), rep(TRUE, 3))
})

test_that("a single floating constraint skips only irrelevant symmetry analysis", {
  ordinary <- "{6S|1,2}Gal(b1-4)Glc"
  symmetric <- "{Neu5Ac(a2-3)|2,4}Gal(a1-?)[Gal(a1-?)]Glc(a1-"
  expected <- .compact_iupac_arrays(c(ordinary, symmetric))
  calls <- 0L
  original <- .compact_bliss_labels
  local_mocked_bindings(.compact_bliss_labels = function(edges, colors) {
    calls <<- calls + 1L
    original(edges, colors)
  })
  expect_identical(.compact_iupac_arrays(ordinary), expected[1])
  expect_identical(calls, 0L)
  expect_identical(.compact_iupac_arrays(symmetric), expected[2])
  expect_identical(calls, 1L)
})
