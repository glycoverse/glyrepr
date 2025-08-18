


test_that("get_mono_type with mixed input works", {
  struct1 <- n_glycan_core(mono_type = "concrete")
  struct2 <- n_glycan_core(mono_type = "generic")
  strucs <- c(struct1, struct2)
  expected <- c("concrete", "generic")
  expect_equal(get_mono_type(strucs), expected)
})


test_that("get mono type of structures", {
  glycan <- n_glycan_core(mono_type = "concrete")
  expect_equal(get_mono_type(glycan), "concrete")
  
  glycan <- n_glycan_core(mono_type = "generic")
  expect_equal(get_mono_type(glycan), "generic")
})


test_that("get mono type of character vector", {
  expect_equal(get_mono_type(c("Gal", "GlcNAc")), c("concrete", "concrete"))
  expect_equal(get_mono_type(c("Hex", "HexNAc")), c("generic", "generic"))
})





test_that("get mono type of composition", {
  comp_generic <- glycan_composition(c(Hex = 4, HexNAc = 1))
  comp_concrete <- glycan_composition(c(Gal = 4, GlcNAc = 1))

  expect_equal(get_mono_type(comp_generic), "generic")
  expect_equal(get_mono_type(comp_concrete), "concrete")
})


# Tests for convert_to_generic function
test_that("convert_to_generic works with character vectors", {
  result <- convert_to_generic(c("Gal", "GlcNAc"))
  expect_equal(result, c("Hex", "HexNAc"))
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
  expect_equal(igraph::V(graph)$mono, c("Hex", "Hex", "Hex", "HexNAc", "HexNAc"))
})


test_that("convert_to_generic with already generic structure returns same", {
  glycan <- n_glycan_core(mono_type = "generic")
  result <- convert_to_generic(glycan)
  expect_identical(result, glycan)
})


test_that("convert_to_generic works with glycan compositions", {
  comp_concrete <- glycan_composition(c(Gal = 2, GlcNAc = 1))
  comp_generic <- convert_to_generic(comp_concrete)

  expect_true(is_glycan_composition(comp_generic))
  data <- vctrs::vec_data(comp_generic)
  expect_equal(vctrs::field(data, "mono_type"), "generic")
  # Check that the conversion aggregated correctly (Gal -> Hex)
  comp_data <- vctrs::field(data, "data")[[1]]
  expect_equal(comp_data, c(Hex = 2, HexNAc = 1))
})


test_that("convert_to_generic with already generic composition returns same", {
  comp_generic <- glycan_composition(c(Hex = 2, HexNAc = 1))
  result <- convert_to_generic(comp_generic)
  expect_identical(result, comp_generic)
})
