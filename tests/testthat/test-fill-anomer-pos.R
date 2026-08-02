test_that("fill_anomer_pos fills missing anomer positions", {
  strucs <- as_glycan_structure(c(
    "Gal(??-?)GalNAc(??-",
    "Neu5Ac(??-?)Gal(??-?)GalNAc(??-"
  ))

  result <- fill_anomer_pos(strucs)

  expect_s3_class(result, "glyrepr_structure")
  expect_equal(
    as.character(result),
    c(
      "Gal(?1-?)GalNAc(?1-",
      "Neu5Ac(?2-?)Gal(?1-?)GalNAc(?1-"
    )
  )
})


test_that("fill_anomer_pos preserves existing anomer annotations", {
  strucs <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    "Neu5Ac(a?-3)Gal(?1-"
  ))

  result <- fill_anomer_pos(strucs)

  expect_equal(
    as.character(result),
    c(
      "Gal(b1-3)GalNAc(a1-",
      "Neu5Ac(a2-3)Gal(?1-"
    )
  )
})

test_that("fill_anomer_pos skips floating normalization for ordinary trees", {
  struc <- as_glycan_structure("Gal(??-?)GalNAc(??-")
  testthat::local_mocked_bindings(
    normalize_floating_parts = function(...) {
      stop("floating normalization should not run")
    }
  )

  result <- fill_anomer_pos(struc)

  expect_identical(
    unname(as.character(result)),
    "Gal(?1-?)GalNAc(?1-"
  )
})


test_that("fill_anomer_pos preserves NA values and names", {
  strucs <- c(
    missing = glycan_structure(NA),
    present = as_glycan_structure("Gal(??-?)GalNAc(??-")
  )

  result <- fill_anomer_pos(strucs)

  expect_equal(names(result), names(strucs))
  expect_true(is.na(result[[1]]))
  expect_equal(as.character(result[[2]]), "Gal(?1-?)GalNAc(?1-")
})


test_that("fill_anomer_pos accepts generic monosaccharides", {
  strucs <- n_glycan_core(linkage = FALSE, mono_type = "generic")

  result <- fill_anomer_pos(strucs)

  expect_s3_class(result, "glyrepr_structure")
  expect_equal(
    as.character(result),
    "Hex(?1-?)[Hex(?1-?)]Hex(?1-?)HexNAc(?1-?)HexNAc(?1-"
  )
})

test_that("fill_anomer_pos fills floating attachment positions", {
  strucs <- as_glycan_structure(c(
    "{Neu5Ac(??-?)|1,2}Gal(??-?)GalNAc(??-",
    "{Gal(??-?)Neu5Ac(??-?)|1,2}GlcNAc(??-?)GalNAc(??-"
  ))

  result <- fill_anomer_pos(strucs)

  expect_identical(
    as.character(result),
    c(
      "{Neu5Ac(?2-?)|1,2}Gal(?1-?)GalNAc(?1-",
      "{Gal(?1-?)Neu5Ac(?2-?)|1,2}GlcNAc(?1-?)GalNAc(?1-"
    )
  )
})

test_that("fill_anomer_pos works with glycan graphs without reordering", {
  structure <- as_glycan_structure(
    "{Neu5Ac(??-?)|1,2}Gal(??-?)GalNAc(??-"
  )
  graph <- get_structure_graphs(structure)
  names_before <- igraph::V(graph)$name
  edges_before <- igraph::as_edgelist(graph, names = FALSE)

  result <- fill_anomer_pos(graph)

  expect_s3_class(result, "igraph")
  expect_identical(igraph::V(result)$name, names_before)
  expect_identical(igraph::as_edgelist(result, names = FALSE), edges_before)
  expect_identical(igraph::E(result)$linkage, "?1-?")
  expect_identical(result$anomer, "?1")
  expect_identical(result$floating_parts[[1]]$linkage, "?2-?")
})

test_that("fill_anomer_pos preserves known floating attachment positions", {
  struc <- as_glycan_structure(
    "{Neu5Ac(a2-3)|1,2}Gal(??-4)GalNAc(??-"
  )

  result <- fill_anomer_pos(struc)

  expect_identical(
    as.character(result),
    "{Neu5Ac(a2-3)|1,2}Gal(?1-4)GalNAc(?1-"
  )
})
