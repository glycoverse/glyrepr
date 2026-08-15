test_that("get mono type of structures", {
  glycan <- n_glycan_core(mono_type = "concrete")
  expect_equal(get_mono_type(glycan), "concrete")

  glycan <- n_glycan_core(mono_type = "generic")
  expect_equal(get_mono_type(glycan), "generic")
})

test_that("get mono type of glycan graphs", {
  concrete <- get_structure_graphs(n_glycan_core(mono_type = "concrete"))
  generic <- get_structure_graphs(n_glycan_core(mono_type = "generic"))

  expect_identical(get_mono_type(concrete), "concrete")
  expect_identical(get_mono_type(generic), "generic")
})

test_that("get mono type of character vector", {
  expect_equal(get_mono_type(c("Gal", "GlcNAc")), c("concrete", "concrete"))
  expect_equal(get_mono_type(c("Hex", "HexNAc")), c("generic", "generic"))
  expect_identical(
    get_mono_type(c(first = "Gal", second = "Hex")),
    c(first = "concrete", second = "generic")
  )
})

test_that("get_mono_type of special monosaccharides", {
  # Some monosaccharides have the same name for both generic and concrete types.
  expect_equal(get_mono_type(c("Mur", "Neu", "Kdn")), rep("concrete", 3))
  expect_equal(get_mono_type(c("gMur", "gNeu", "gKdn")), rep("generic", 3))
})

test_that("get_mono_type of unknown monosaccharides", {
  expect_error(
    get_mono_type(c("Unknown", "Unknown2")),
    "Unknown monosaccharide"
  )
})

test_that("get mono type of composition", {
  comp_generic <- glycan_composition(
    c(Hex = 4, HexNAc = 1),
    c(Hex = 1, HexNAc = 1)
  )
  comp_concrete <- glycan_composition(c(Gal = 4, GlcNAc = 1))

  expect_equal(get_mono_type(comp_generic), c("generic", "generic"))
  expect_equal(get_mono_type(comp_concrete), "concrete")
})

test_that("get_mono_type reports mixed elements and preserves names", {
  structures <- as_glycan_structure(c(
    concrete = "Gal(b1-3)GalNAc(a1-",
    generic = "Hex(b1-3)HexNAc(a1-",
    mixed = "Hex(b1-3)GalNAc(a1-",
    missing = NA_character_
  ))
  compositions <- glycan_composition(
    concrete = c(Gal = 1, GalNAc = 1),
    generic = c(Hex = 1, HexNAc = 1),
    mixed = c(Hex = 1, GalNAc = 1),
    missing = NA
  )

  expected <- c(
    concrete = "concrete",
    generic = "generic",
    mixed = "mixed",
    missing = NA_character_
  )
  expect_identical(get_mono_type(structures), expected)
  expect_identical(get_mono_type(compositions), unname(expected))
  expect_identical(
    get_mono_type(get_structure_graphs(structures[[3]])),
    "mixed"
  )
})

test_that("get_mono_type of composition with substituents", {
  # Concrete composition with sulfate substituent
  comp_sulfate <- glycan_composition(c(Gal = 2, GlcNAc = 1, S = 1))
  expect_equal(get_mono_type(comp_sulfate), "concrete")

  # Generic composition with methyl substituent
  comp_methyl <- glycan_composition(c(Hex = 2, HexNAc = 1, Me = 1))
  expect_equal(get_mono_type(comp_methyl), "generic")

  # Composition with multiple different substituents
  comp_multi_sub <- glycan_composition(c(Gal = 3, GlcNAc = 2, S = 1, Ac = 1))
  expect_equal(get_mono_type(comp_multi_sub), "concrete")
})

test_that("get mono type of empty composition", {
  comp_empty <- glycan_composition()
  expect_equal(get_mono_type(comp_empty), character())
})

# Tests for convert_to_generic function
test_that("convert_to_generic works with character vectors", {
  result <- convert_to_generic(c("Gal", "GlcNAc"))
  expect_equal(result, c("Hex", "HexNAc"))
})

test_that("convert_to_generic converts concrete elements in mixed character vectors", {
  result <- convert_to_generic(c("Gal", "GalNAc", "Hex", "HexNAc"))
  expect_equal(result, c("Hex", "HexNAc", "Hex", "HexNAc"))
})

test_that("convert_to_generic with already generic characters returns same", {
  input <- c("Hex", "HexNAc")
  result <- convert_to_generic(input)
  expect_identical(result, input)
})

test_that("convert_to_generic works with glycan structures", {
  glycan <- n_glycan_core(mono_type = "concrete")
  glycan_generic <- convert_to_generic(glycan)

  expect_true(is_glycan_structure(glycan_generic))
  graph <- get_structure_graphs(glycan_generic, return_list = FALSE)
  expect_equal(
    igraph::V(graph)$mono,
    c("Hex", "Hex", "Hex", "HexNAc", "HexNAc")
  )
})

test_that("convert_to_generic converts mixed structures", {
  glycan <- as_glycan_structure("Hex(b1-3)GalNAc(a1-")
  result <- convert_to_generic(glycan)

  expect_identical(as.character(result), "Hex(b1-3)HexNAc(a1-")
  expect_identical(get_mono_type(result), "generic")
})

test_that("convert_to_generic handles mixed structure vectors", {
  structures <- as_glycan_structure(c(
    concrete = "Gal(b1-3)GalNAc(a1-",
    generic = "Hex(b1-3)HexNAc(a1-",
    mixed = "Hex(b1-3)GalNAc(a1-",
    missing = NA_character_
  ))

  result <- convert_to_generic(structures)

  expect_identical(
    get_mono_type(result),
    c(
      concrete = "generic",
      generic = "generic",
      mixed = "generic",
      missing = NA_character_
    )
  )
  expect_identical(names(result), names(structures))
})

test_that("convert_to_generic works with glycan graphs without reordering", {
  structure <- as_glycan_structure(
    "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
  )
  graph <- get_structure_graphs(structure)
  names_before <- igraph::V(graph)$name
  edges_before <- igraph::as_edgelist(graph, names = FALSE)

  result <- convert_to_generic(graph)

  expect_s3_class(result, "igraph")
  expect_identical(igraph::V(result)$name, names_before)
  expect_identical(igraph::as_edgelist(result, names = FALSE), edges_before)
  expect_identical(igraph::V(result)$mono, c("NeuAc", "Hex", "HexNAc"))
  expect_identical(result$floating_parts, graph$floating_parts)
})

test_that("convert_to_generic with already generic structure returns same", {
  glycan <- n_glycan_core(mono_type = "generic")
  result <- convert_to_generic(glycan)
  expect_identical(result, glycan)
})

test_that("convert_to_generic preserves empty structure vectors", {
  glycan <- as_glycan_structure(character())

  expect_identical(convert_to_generic(glycan), glycan)
})

test_that("convert_to_generic works with glycan compositions", {
  comp_concrete <- glycan_composition(c(Gal = 2, GlcNAc = 1))
  comp_generic <- convert_to_generic(comp_concrete)

  expect_true(is_glycan_composition(comp_generic))
  expect_equal(as.character(comp_generic), "Hex(2)HexNAc(1)")
})

test_that("convert_to_generic converts and aggregates mixed compositions", {
  comp <- glycan_composition(c(Hex = 2, Gal = 1, Glc = 3, GalNAc = 1))
  result <- convert_to_generic(comp)

  expect_identical(as.character(result), "Hex(6)HexNAc(1)")
  expect_identical(get_mono_type(result), "generic")
})

test_that("convert_to_generic handles mixed composition vectors", {
  comps <- glycan_composition(
    c(Gal = 1),
    c(Hex = 1),
    c(Hex = 1, Gal = 2),
    NA
  )

  result <- convert_to_generic(comps)

  expect_identical(
    get_mono_type(result),
    c("generic", "generic", "generic", NA_character_)
  )
  expect_identical(as.character(result[1:3]), c("Hex(1)", "Hex(1)", "Hex(3)"))
  expect_identical(is.na(result), c(FALSE, FALSE, FALSE, TRUE))
})

test_that("convert_to_generic works with glycan compositions with substituents", {
  comp <- glycan_composition(c(Gal = 2, GlcNAc = 1, Me = 1))
  result <- convert_to_generic(comp)
  expected <- glycan_composition(c(Hex = 2, HexNAc = 1, Me = 1))
  expect_equal(result, expected)
})

test_that("convert_to_generic with already generic composition returns same", {
  comp_generic <- glycan_composition(c(Hex = 2, HexNAc = 1))
  result <- convert_to_generic(comp_generic)
  expect_identical(result, comp_generic)
})

test_that("convert_to_generic with empty composition returns empty composition", {
  comp_empty <- glycan_composition()
  result <- convert_to_generic(comp_empty)
  expect_equal(result, comp_empty)
})

test_that("convert_to_generic preserves NA in compositions", {
  # NA as second element
  comps1 <- glycan_composition(c(Gal = 1), NA)
  result1 <- convert_to_generic(comps1)
  expect_equal(length(result1), 2)
  expect_equal(as.character(result1[1]), "Hex(1)")
  expect_true(is.na(result1[2]))

  # NA as first element
  comps2 <- glycan_composition(NA, c(Gal = 1))
  result2 <- convert_to_generic(comps2)
  expect_equal(length(result2), 2)
  expect_true(is.na(result2[1]))
  expect_equal(as.character(result2[2]), "Hex(1)")

  # NA in the middle
  comps3 <- glycan_composition(c(Gal = 1), NA, c(GlcNAc = 2))
  result3 <- convert_to_generic(comps3)
  expect_equal(length(result3), 3)
  expect_equal(as.character(result3[1]), "Hex(1)")
  expect_true(is.na(result3[2]))
  expect_equal(as.character(result3[3]), "HexNAc(2)")
})

test_that("convert_to_generic handles NA in glyrepr_structure", {
  # Create structure vector with NA
  # Use valid IUPAC-condensed strings from existing structures
  glycans <- as_glycan_structure(c(
    "Gal(b1-3)GalNAc(a1-",
    NA,
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))

  result <- convert_to_generic(glycans)
  expect_equal(length(result), 3)
  expect_false(is.na(result[1]))
  expect_true(is.na(result[2]))
  expect_false(is.na(result[3]))

  # Verify non-NA elements are converted
  expect_true(all(grepl("Hex", as.character(result[!is.na(result)]))))
})

test_that("convert_to_generic with all NA structures", {
  # All NA elements
  glycans <- as_glycan_structure(c(NA, NA))
  result <- convert_to_generic(glycans)
  expect_equal(length(result), 2)
  expect_true(all(is.na(result)))
})

test_that("get_mono_type aligns NA compositions", {
  concrete <- glycan_composition(NA, c(Gal = 1))
  generic <- glycan_composition(NA, c(Hex = 1))
  all_na <- glycan_composition(NA, NA)

  expect_identical(get_mono_type(concrete), c(NA, "concrete"))
  expect_identical(get_mono_type(generic), c(NA, "generic"))
  expect_identical(get_mono_type(all_na), c(NA_character_, NA_character_))
})
