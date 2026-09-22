test_that("singleton localization prunes occupied candidate domains", {
  inputs <- c(
    "{Neu5Ac(a1-3)|4}{Fuc(a1-3)|3,5,4}Gal6S(b1-4)[Gal(b1-?)]Glc",
    paste0(
      "{3/6Ac|1,4}{Gal6S(b1-3/6)}{Neu5Ac(a1-3)|4}",
      "{Gal(b1-3/6)|6,1}Gal6S(b1-4)[Gal(b1-?)]Glc"
    ),
    "{Fuc(a1-3)|3}{Gal(b1-3)|3,4}Glc(b1-4)Man",
    "{3Ac|2}{Gal(b1-3)|2,3}Glc(b1-4)Man",
    "{3Ac|2,3}{Fuc(a1-3)|2}Glc(b1-4)Man",
    "{3Ac|2}{Fuc(a1-3)}Glc(b1-4)Man(b1-3)Gal",
    "{3Ac}{3S|1}Glc(b1-4)Man"
  )

  for (input in inputs) {
    raw <- .parse_iupac_condensed_single(input)
    expect_no_error(validate_glycan_graph(raw))
    reference <- canonicalize_floating_graph(raw)
    expect_no_error(validate_glycan_graph(reference))

    expect_identical(.compact_iupac_arrays(input)[[1]]$status, "ok")
    structure <- as_glycan_structure(input)
    graph <- get_structure_graphs(structure, return_list = FALSE)
    expect_no_error(validate_glycan_graph(graph))
    raw_variants <- enumerate_floating_graph_localizations(raw)$graph
    normalized_variants <- enumerate_floating_graph_localizations(graph)$graph
    expect_setequal(
      vapply(raw_variants, graph_to_iupac, character(1)),
      vapply(normalized_variants, graph_to_iupac, character(1))
    )
    expect_identical(as.character(structure), graph_to_iupac(reference))
    expect_identical(
      as.character(as_glycan_structure(as.character(structure))),
      as.character(structure)
    )
    expect_identical(
      as.character(as_glycan_structure(raw)),
      as.character(structure)
    )
  }
})


test_that("unrestricted domains resolve when localization leaves one parent", {
  inputs <- c(
    "{3Ac|2}{Fuc(a1-3)}Glc(b1-4)Man(b1-3)Gal",
    "{3Ac}{3S|1}Glc(b1-4)Man"
  )
  expected <- as_glycan_structure(c(
    "Glc3Ac(b1-4)[Fuc(a1-3)]Man(b1-3)Gal",
    "Glc3S(b1-4)Man3Ac"
  ))
  expect_identical(
    as.character(as_glycan_structure(inputs)),
    as.character(expected)
  )
  for (i in seq_along(inputs)) {
    graph <- canonicalize_floating_graph(.parse_iupac_condensed_single(inputs[
      i
    ]))
    expect_identical(has_floating_metadata(graph), FALSE)
    expect_identical(graph_to_iupac(graph), as.character(expected[i]))
  }
})
