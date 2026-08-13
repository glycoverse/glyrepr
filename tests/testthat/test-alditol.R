test_that("alditol IUPAC syntax round-trips", {
  iupacs <- c(
    single = "Glc-ol(a1-",
    branched = "Gal(a1-3)[Fuc(a1-4)]GlcNAc-ol(b1-",
    substituted = "Gal(b1-4)GlcNAc6S-ol(a1-",
    inferred = "Glc-ol",
    floating = "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc-ol(a1-"
  )

  glycans <- as_glycan_structure(iupacs)

  expect_identical(
    unname(structure_to_iupac(glycans)),
    unname(c(iupacs[1:3], "Glc-ol(?1-", iupacs[5]))
  )
  expect_identical(unname(get_alditol(glycans)), rep(TRUE, 5))
  expect_identical(
    unname(purrr::map_lgl(as.list(glycans), ~ .x$alditol)),
    rep(TRUE, 5)
  )
  expect_identical(
    structure_to_iupac(get_structure_graphs(glycans[[1]])),
    "Glc-ol(a1-"
  )
})


test_that("ordinary glycans default to non-alditols", {
  graph <- get_structure_graphs(o_glycan_core_1())
  legacy <- igraph::delete_graph_attr(graph, "alditol")

  expect_identical(validate_glycan_graph(legacy), legacy)
  expect_identical(get_alditol(legacy), FALSE)

  canonical <- canonicalize_glycan_graph(legacy)
  expect_identical(igraph::graph_attr(canonical, "alditol"), FALSE)
  expect_identical(graph_to_iupac(canonical), "Gal(b1-3)GalNAc(a1-")
})


test_that("alditol graph attributes are validated", {
  graph <- get_structure_graphs(o_glycan_core_1())

  graph$alditol <- TRUE
  expect_identical(validate_glycan_graph(graph), graph)

  graph$alditol <- FALSE
  expect_identical(validate_glycan_graph(graph), graph)

  graph$alditol <- "TRUE"
  expect_snapshot(validate_glycan_graph(graph), error = TRUE)

  graph$alditol <- c(TRUE, FALSE)
  expect_snapshot(validate_glycan_graph(graph), error = TRUE)

  graph$alditol <- NA
  expect_snapshot(validate_glycan_graph(graph), error = TRUE)
})


test_that("alditol status distinguishes canonical structure keys", {
  glycans <- as_glycan_structure(c(
    reduced = "Glc-ol(a1-",
    ordinary = "Glc(a1-",
    duplicate = "Glc-ol(a1-",
    missing = NA
  ))

  expect_identical(
    unname(get_alditol(glycans)),
    c(TRUE, FALSE, TRUE, NA)
  )
  expect_identical(names(get_alditol(glycans)), names(glycans))
  expect_length(attr(glycans, "graphs"), 2)
  expect_identical(
    unname(structure_to_iupac(glycans)),
    c("Glc-ol(a1-", "Glc(a1-", "Glc-ol(a1-", NA)
  )
})


test_that("alditol status survives structure transformations", {
  glycan <- as_glycan_structure("Gal6S(b1-4)GlcNAc-ol(a1-")
  transformed <- list(
    convert_to_generic(glycan),
    fill_anomer_pos(glycan),
    remove_linkages(glycan),
    remove_substituents(glycan),
    smap_structure(glycan, identity)
  )

  expect_identical(
    purrr::map_lgl(transformed, ~ unname(get_alditol(.x))),
    rep(TRUE, length(transformed))
  )

  normalized <- smap_structure(
    as_glycan_structure("Glc(a1-"),
    ~ igraph::delete_graph_attr(.x, "alditol")
  )
  expect_identical(
    igraph::graph_attr(get_structure_graphs(normalized), "alditol"),
    FALSE
  )

  floating <- as_glycan_structure(
    "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc-ol(a1-"
  )
  variants <- enumerate_floating_localizations(floating)
  expect_identical(
    unname(get_alditol(variants$structure)),
    rep(TRUE, nrow(variants))
  )
  expect_identical(
    purrr::map_lgl(
      enumerate_floating_graph_localizations(
        get_structure_graphs(floating)
      )$graph,
      get_alditol
    ),
    rep(TRUE, 2)
  )
})


test_that("alditol markers are restricted to the main reducing end", {
  expect_snapshot(
    as_glycan_structure("Gal-ol(b1-4)Glc(a1-"),
    error = TRUE
  )
  expect_snapshot(
    as_glycan_structure("Glc-ol-ol(a1-"),
    error = TRUE
  )
  expect_snapshot(
    as_glycan_structure("{Gal-ol(b1-4)|1}Glc(a1-"),
    error = TRUE
  )
})


test_that("alditol formatting colors the reducing-end monosaccharide", {
  glycan <- as_glycan_structure("Gal(b1-4)GlcNAc-ol(a1-")

  expect_identical(
    cli::ansi_strip(format_glycan_structure_subset(glycan, 1)),
    "Gal(b1-4)GlcNAc-ol(a1-"
  )
})
