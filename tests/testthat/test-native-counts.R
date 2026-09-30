test_that("native graph counts preserve reference compositions and query semantics", {
  graphs <- lapply(available_monosaccharides(), function(mono) {
    g <- igraph::make_empty_graph(2, directed = TRUE)
    g <- igraph::set_vertex_attr(g, "mono", value = c(mono, "Gal"))
    igraph::set_vertex_attr(g, "sub", value = c("3Me,6S", ""))
  })
  reference <- lapply(graphs, .graph_to_composition_reference)
  expect_identical(.native_graph_counts(graphs), reference)
  comps <- new_glycan_composition(reference)
  for (query in list(
    NULL,
    "Hex",
    "HexNAc",
    "dHex",
    "Gal",
    "D-Fuc",
    "Galf",
    "Kdn",
    "S",
    "Me"
  )) {
    for (subs in c(FALSE, TRUE)) {
      expect_identical(
        .native_graph_counts(graphs, query, subs, TRUE),
        count_mono(comps, query, subs)
      )
    }
  }
})

test_that("structure counts preserve duplicates, missing values and names", {
  x <- c(n_glycan_core(), o_glycan_core_1(), n_glycan_core(), NA)
  names(x) <- letters[seq_along(x)]
  reference <- new_glycan_composition(smap(x, .graph_to_composition_reference))
  expect_identical(as_glycan_composition(x), reference)
  for (query in list(NULL, "Hex", "Gal", "Me")) {
    expect_identical(count_mono(x, query), count_mono(reference, query))
  }
  expect_identical(count_mono(glycan_structure()), integer())
  expect_identical(count_mono(x[4]), c(d = NA_integer_))
})

test_that("unsupported substituent representations use the reference", {
  g <- igraph::make_empty_graph(1, directed = TRUE)
  g <- igraph::set_vertex_attr(g, "mono", value = "Gal")
  g <- igraph::set_vertex_attr(g, "sub", value = "3Unknown")
  expect_null(.native_graph_counts(list(g)))
  expect_identical(graph_to_composition(g), .graph_to_composition_reference(g))
})

test_that("floating substituents and parts contribute once", {
  x <- as_glycan_structure(c(
    "{6S}Gal(a1-3)Gal(a1-",
    "{?S}Gal(a1-3)Glc(a1-3)Man(a1-"
  ))
  reference <- new_glycan_composition(smap(x, .graph_to_composition_reference))
  expect_identical(as_glycan_composition(x), reference)
  expect_identical(count_mono(x, "S"), c(1L, 1L))
  expect_identical(count_mono(x, include_subs = TRUE), c(3L, 4L))
})
