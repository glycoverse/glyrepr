test_that("get_structure_level is element-wise and ignores residue type", {
  glycans <- as_glycan_structure(c(
    concrete_intact = "Gal(b1-3)GalNAc(a1-",
    generic_intact = "Hex(b1-3)HexNAc(a1-",
    mixed_intact = "Hex(b1-3)GalNAc(a1-",
    concrete_partial = "Gal(b1-?)GalNAc(a1-",
    generic_partial = "Hex(b1-?)HexNAc(a1-",
    mixed_partial = "Hex(b1-?)GalNAc(a1-",
    concrete_topological = "Gal(??-?)GalNAc(??-",
    generic_topological = "Hex(??-?)HexNAc(??-",
    mixed_topological = "Hex(??-?)GalNAc(??-",
    missing = NA_character_
  ))

  expect_identical(
    get_structure_level(glycans),
    c(
      concrete_intact = "intact",
      generic_intact = "intact",
      mixed_intact = "intact",
      concrete_partial = "partial",
      generic_partial = "partial",
      mixed_partial = "partial",
      concrete_topological = "topological",
      generic_topological = "topological",
      mixed_topological = "topological",
      missing = NA_character_
    )
  )
})

test_that("get_structure_level works with glycan graphs", {
  structures <- as_glycan_structure(c(
    "Hex(b1-3)GalNAc(a1-",
    "Hex(b1-?)GalNAc(a1-",
    "Hex(??-?)GalNAc(??-"
  ))
  graphs <- get_structure_graphs(structures)

  expect_identical(
    vapply(graphs, get_structure_level, character(1)),
    c("intact", "partial", "topological")
  )
})

test_that("ambiguous linkage positions are partial", {
  glycan <- as_glycan_structure("Neu5Ac(a2-3/6)Gal(b1-4)GlcNAc(b1-")
  expect_identical(get_structure_level(glycan), "partial")
})

test_that("reducing-end information contributes to structure level", {
  glycans <- as_glycan_structure(c(
    "Gal(??-?)GalNAc(a1-",
    "Gal(??-?)GalNAc(?1-",
    "Gal(??-?)GalNAc(??-"
  ))

  expect_identical(
    get_structure_level(glycans),
    c("partial", "partial", "topological")
  )
})

test_that("structure level includes floating attachment linkages", {
  glycans <- as_glycan_structure(c(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-",
    "{NeuAc(a2-?)|2,3}Hex(b1-3)GalNAc(a1-",
    "{NeuAc(??-?)|2,3}Hex(??-?)GalNAc(??-"
  ))

  expect_identical(
    get_structure_level(glycans),
    c("intact", "partial", "topological")
  )
})

test_that("floating parent ambiguity does not change structure level", {
  intact <- as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )
  topological <- as_glycan_structure(
    "{NeuAc(??-?)|2,3}Hex(??-?)GalNAc(??-"
  )

  expect_identical(get_structure_level(intact), "intact")
  expect_identical(get_structure_level(topological), "topological")
})

test_that("floating substituents do not change structure level", {
  glycans <- as_glycan_structure(c(
    "{6S|1,2}Gal(b1-3)GalNAc(a1-",
    "{?S|1,2}Hex(??-?)GalNAc(??-"
  ))

  expect_identical(
    get_structure_level(glycans),
    c("intact", "topological")
  )
})

test_that("get_structure_level preserves missingness and names", {
  all_na <- as_glycan_structure(c(
    first = NA_character_,
    second = NA_character_
  ))
  empty <- as_glycan_structure(character())

  expect_identical(
    get_structure_level(all_na),
    c(first = NA_character_, second = NA_character_)
  )
  expect_identical(get_structure_level(empty), character())
})

test_that("remove_linkages produces topological structures", {
  glycans <- as_glycan_structure(c(
    concrete = "Gal(b1-3)GalNAc(a1-",
    generic = "Hex(b1-3)HexNAc(a1-",
    mixed = "Hex(b1-3)GalNAc(a1-",
    missing = NA_character_
  ))

  result <- remove_linkages(glycans)

  expect_identical(
    get_structure_level(result),
    c(
      concrete = "topological",
      generic = "topological",
      mixed = "topological",
      missing = NA_character_
    )
  )
  expect_identical(get_mono_type(result), get_mono_type(glycans))
})

test_that("remove_linkages preserves graph order and floating metadata", {
  structure <- as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )
  graph <- get_structure_graphs(structure)
  names_before <- igraph::V(graph)$name
  edges_before <- igraph::as_edgelist(graph, names = FALSE)

  result <- remove_linkages(graph)

  expect_s3_class(result, "igraph")
  expect_identical(igraph::V(result)$name, names_before)
  expect_identical(igraph::as_edgelist(result, names = FALSE), edges_before)
  expect_identical(igraph::E(result)$linkage, "??-?")
  expect_identical(result$anomer, "??")
  expect_identical(result$floating_parts[[1]]$linkage, "??-?")
  expect_identical(get_structure_level(result), "topological")
})
