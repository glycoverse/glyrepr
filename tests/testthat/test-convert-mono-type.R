test_that("convert from concrete to generic for NE graph", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "concrete")

  glycan_generic <- convert_glycan_mono_type(glycan, to = "generic")

  expect_equal(igraph::V(glycan_generic)$mono, c("HexNAc", "HexNAc", "Hex", "Hex", "Hex"))
})


test_that("convert from concrete to generic for DN graph", {
  glycan <- n_glycan_core(mode = "dn", mono_type = "concrete")

  glycan_generic <- convert_glycan_mono_type(glycan, to = "generic")

  monos <- igraph::V(glycan_generic)$mono
  monos <- monos[!is.na(monos)]
  expect_equal(monos, c("HexNAc", "HexNAc", "Hex", "Hex", "Hex"))
})


test_that("convert from generic to simple for NE graph", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")

  glycan_simple <- convert_glycan_mono_type(glycan, to = "simple")

  expect_equal(igraph::V(glycan_simple)$mono, c("N", "N", "H", "H", "H"))
})


test_that("convert from generic to simple for DN graph", {
  glycan <- n_glycan_core(mode = "dn", mono_type = "generic")

  glycan_simple <- convert_glycan_mono_type(glycan, to = "simple")

  monos <- igraph::V(glycan_simple)$mono
  monos <- monos[!is.na(monos)]
  expect_equal(monos, c("N", "N", "H", "H", "H"))
})


test_that("convert from concrete to simple for NE graph", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "concrete")

  glycan_simple <- convert_glycan_mono_type(glycan, to = "simple")

  expect_equal(igraph::V(glycan_simple)$mono, c("N", "N", "H", "H", "H"))
})


test_that("convert from concrete to simple for DN graph", {
  glycan <- n_glycan_core(mode = "dn", mono_type = "concrete")

  glycan_simple <- convert_glycan_mono_type(glycan, to = "simple")

  monos <- igraph::V(glycan_simple)$mono
  monos <- monos[!is.na(monos)]
  expect_equal(monos, c("N", "N", "H", "H", "H"))
})


test_that("converting from simple to generic fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "generic"), error = TRUE)
})


test_that("converting from simple to concrete fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("converting from generic to concrete fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("converting from generic to generic fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "generic"), error = TRUE)
})


test_that("converting from simple to simple fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "simple"), error = TRUE)
})


test_that("converting from concrete to concrete fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "concrete")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "concrete"), error = TRUE)
})
