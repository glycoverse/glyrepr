test_that("convert from concrete to generic", {
  glycan <- n_glycan_core(mono_type = "concrete")

  glycan_generic <- convert_mono_type(glycan, to = "generic")
  
  # Extract the igraph from the vectorized structure
  graph <- get_structure_graphs(glycan_generic, 1)
  expect_equal(igraph::V(graph)$mono, c("HexNAc", "HexNAc", "Hex", "Hex", "Hex"))
})


test_that("convert from generic to simple", {
  glycan <- n_glycan_core(mono_type = "generic")

  glycan_simple <- convert_mono_type(glycan, to = "simple")

  # Extract the igraph from the vectorized structure
  graph <- get_structure_graphs(glycan_simple, 1)
  expect_equal(igraph::V(graph)$mono, c("N", "N", "H", "H", "H"))
})


test_that("convert from concrete to simple", {
  glycan <- n_glycan_core(mono_type = "concrete")

  glycan_simple <- convert_mono_type(glycan, to = "simple")

  # Extract the igraph from the vectorized structure
  graph <- get_structure_graphs(glycan_simple, 1)
  expect_equal(igraph::V(graph)$mono, c("N", "N", "H", "H", "H"))
})


test_that("converting glycan from simple to generic fails", {
  glycan <- n_glycan_core(mono_type = "simple")
  expect_snapshot(convert_mono_type(glycan, to = "generic"), error = TRUE)
})


test_that("converting glycan from simple to concrete fails", {
  glycan <- n_glycan_core(mono_type = "simple")
  expect_snapshot(convert_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("converting glycan from generic to concrete fails", {
  glycan <- n_glycan_core(mono_type = "generic")
  expect_snapshot(convert_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("converting glycan from generic to generic returns same object", {
  glycan <- n_glycan_core(mono_type = "generic")
  result <- convert_mono_type(glycan, to = "generic")
  expect_identical(result, glycan)
})


test_that("converting glycan from simple to simple returns same object", {
  glycan <- n_glycan_core(mono_type = "simple")
  result <- convert_mono_type(glycan, to = "simple")
  expect_identical(result, glycan)
})


test_that("converting glycan mono types with NA produced", {
  glycan_vec <- o_glycan_core_1()
  glycan <- get_structure_graphs(glycan_vec, 1)  # Extract the igraph
  igraph::V(glycan)$mono[[1]] <- "Pse"  # cannot be converted to generic
  
  # Create a new vectorized structure with the modified graph
  modified_glycan_vec <- glycan_structure(glycan)
  expect_snapshot(convert_mono_type(modified_glycan_vec, to = "generic"), error = TRUE)
})


test_that("convert mono types", {
  before <- c("Man", "Hex", "GlcNAc", "dHex", "Neu5Ac", "Neu5Gc", "Sia")
  to <- "simple"
  after <- c("H", "H", "N", "F", "A", "G", "S")
  expect_equal(convert_mono_type(before, to), after)
})


test_that("convert mono types returns same for already target type", {
  before <- c("H", "N", "F")
  to <- "simple"
  result <- convert_mono_type(before, to)
  expect_equal(result, before)
})


test_that("convert mono types with mixed types", {
  before <- c("H", "Hex", "N")  # mixed simple and generic
  to <- "simple"
  after <- c("H", "H", "N")
  expect_equal(convert_mono_type(before, to), after)
})


test_that("convert mono types with bad directions", {
  before <- c("H", "Hex")
  to <- "concrete"
  expect_snapshot(convert_mono_type(before, to), error = TRUE)
})


test_that("converting mono types with NA", {
  # This should throw a warning
  expect_snapshot(convert_mono_type(c("Pse", "Fuc"), "generic"))
})


patrick::with_parameters_test_that("get glycan mono type works", {
    expect_equal(get_mono_type(glycan), type)
  },
  glycan = list(
    n_glycan_core(mono_type = "concrete"),
    n_glycan_core(mono_type = "generic"),
    n_glycan_core(mono_type = "simple")
  ),
  type = c("concrete", "generic", "simple")
)


test_that("get mono types vectorized", {
  mono <- c("Gal", "Hex", "H")
  expected <- c("concrete", "generic", "simple")
  expect_equal(get_mono_type(mono), expected)
})


test_that("get mono types fails for unknown monos", {
  mono <- c("bad1", "bad2", "bad3")
  expect_snapshot(get_mono_type(mono), error = TRUE)
})


# Test glyrepr_composition conversion
test_that("convert composition mono types", {
  comp_concrete <- glycan_composition(c(Glc = 2, GalNAc = 1))
  comp_generic <- convert_mono_type(comp_concrete, to = "generic")
  
  # Check that it converts correctly
  data <- vctrs::vec_data(comp_generic)
  comp_data <- vctrs::field(data, "data")[[1]]
  expect_equal(names(comp_data), c("Hex", "HexNAc"))
  expect_equal(as.integer(comp_data), c(2L, 1L))
})


test_that("convert composition to simple", {
  comp_generic <- glycan_composition(c(Hex = 3, HexNAc = 2, dHex = 1))
  comp_simple <- convert_mono_type(comp_generic, to = "simple")
  
  # Check that it converts correctly
  data <- vctrs::vec_data(comp_simple)
  comp_data <- vctrs::field(data, "data")[[1]]
  expect_equal(names(comp_data), c("H", "N", "F"))
  expect_equal(as.integer(comp_data), c(3L, 2L, 1L))
})


test_that("convert composition returns same for target type", {
  comp_simple <- glycan_composition(c(H = 2, N = 1))
  result <- convert_mono_type(comp_simple, to = "simple")
  expect_identical(result, comp_simple)
})


test_that("convert composition with aggregation", {
  # Create a composition where conversion would result in duplicate mono types
  comp_concrete <- glycan_composition(c(Glc = 2, Man = 1, Gal = 1))
  comp_simple <- convert_mono_type(comp_concrete, to = "simple")
  
  # All should become "H" and aggregate to 4
  data <- vctrs::vec_data(comp_simple)
  comp_data <- vctrs::field(data, "data")[[1]]
  expect_equal(names(comp_data), "H")
  expect_equal(as.integer(comp_data), 4L)
})


# Test get_mono_type for compositions
test_that("get mono type for composition", {
  comp_concrete <- glycan_composition(c(Glc = 2, GalNAc = 1))
  comp_generic <- glycan_composition(c(Hex = 3, HexNAc = 2))
  comp_simple <- glycan_composition(c(H = 4, N = 1))
  
  expect_equal(get_mono_type(comp_concrete), "concrete")
  expect_equal(get_mono_type(comp_generic), "generic")
  expect_equal(get_mono_type(comp_simple), "simple")
})


test_that("get mono type for multiple compositions", {
  comp1 <- glycan_composition(c(Glc = 2, GalNAc = 1), c(Man = 1, GlcNAc = 2))
  comp2 <- glycan_composition(c(Hex = 3, HexNAc = 2), c(Hex = 1, dHex = 1))
  
  expect_equal(get_mono_type(comp1), c("concrete", "concrete"))
  expect_equal(get_mono_type(comp2), c("generic", "generic"))
})
