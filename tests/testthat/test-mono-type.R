test_that("convert structure from concrete to generic", {
  glycan <- n_glycan_core(mono_type = "concrete")
  glycan_generic <- convert_mono_type(glycan, to = "generic")
  
  expect_true(is_glycan_structure(glycan_generic))
  graph <- get_structure_graphs(glycan_generic, 1)
  expect_equal(igraph::V(graph)$mono, c("HexNAc", "HexNAc", "Hex", "Hex", "Hex"))
})


test_that("converting glycan from generic to concrete fails", {
  glycan <- n_glycan_core(mono_type = "generic")
  expect_error(convert_mono_type(glycan, to = "concrete"))
})


test_that("converting character vectors from concrete to generic", {
  result <- convert_mono_type(c("Gal", "GlcNAc"), to = "generic")
  expect_equal(result, c("Hex", "HexNAc"))
})


test_that("converting character vectors from generic to concrete fails", {
  result <- tryCatch(
    convert_mono_type(c("Hex", "HexNAc"), to = "concrete"),
    error = function(e) e$message
  )
  expect_match(result, "cannot be converted")
})


test_that("converting glycan from concrete to concrete returns same object", {
  glycan <- n_glycan_core(mono_type = "concrete")
  result <- convert_mono_type(glycan, to = "concrete")
  expect_identical(glycan, result)
})


test_that("converting glycan from generic to generic returns same object", {
  glycan <- n_glycan_core(mono_type = "generic")
  result <- convert_mono_type(glycan, to = "generic")
  expect_identical(glycan, result)
})


test_that("convert character: concrete to generic", {
  before <- c("Gal", "GlcNAc", "Fuc")
  to <- "generic"
  result <- convert_mono_type(before, to)
  expect_equal(result, c("Hex", "HexNAc", "dHex"))
})


test_that("convert character: generic to concrete (should fail)", {
  before <- c("Hex", "HexNAc", "dHex")
  to <- "concrete"
  expect_error(convert_mono_type(before, to))
})


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


test_that("convert composition to generic", {
  comp_concrete <- glycan_composition(c(Gal = 2, GlcNAc = 1))
  comp_generic <- convert_mono_type(comp_concrete, to = "generic")
  
  expect_true(is_glycan_composition(comp_generic))
  data <- vctrs::vec_data(comp_generic)
  expect_equal(vctrs::field(data, "mono_type"), "generic")
  # Check that the conversion aggregated correctly (Gal -> Hex)
  comp_data <- vctrs::field(data, "data")[[1]]
  expect_equal(comp_data, c(Hex = 2, HexNAc = 1))
})


test_that("convert composition to same type returns same object", {
  comp_generic <- glycan_composition(c(Hex = 2, HexNAc = 1))
  result <- convert_mono_type(comp_generic, to = "generic")
  expect_identical(result, comp_generic)
})


test_that("convert composition from generic to concrete fails", {
  comp_generic <- glycan_composition(c(Hex = 2, HexNAc = 1))
  expect_error(convert_mono_type(comp_generic, to = "concrete"))
})


test_that("get mono type of composition", {
  comp_generic <- glycan_composition(c(Hex = 4, HexNAc = 1))
  comp_concrete <- glycan_composition(c(Gal = 4, GlcNAc = 1))
  
  expect_equal(get_mono_type(comp_generic), "generic")
  expect_equal(get_mono_type(comp_concrete), "concrete")
})
