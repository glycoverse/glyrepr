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


test_that("converting glycan from simple to generic fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "generic"), error = TRUE)
})


test_that("converting glycan from simple to concrete fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("converting glycan from generic to concrete fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("converting glycan from generic to generic fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "generic"), error = TRUE)
})


test_that("converting glycan from simple to simple fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "simple"), error = TRUE)
})


test_that("converting glycan from concrete to concrete fails", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "concrete")
  expect_snapshot(convert_glycan_mono_type(glycan, to = "concrete"), error = TRUE)
})


test_that("convert mono types", {
  before = c("Man", "Hex", "GlcNAc")
  to = "simple"
  after = c("H", "H", "N")
  expect_equal(convert_mono_type(before, to), after)
})


test_that("convert mono types fails for monos already in simple form", {
  before = c("H", "H", "N", "Hex")
  to = "simple"
  expect_snapshot(convert_mono_type(before, to), error = TRUE)
})


test_that("convert mono types with bad directions", {
  before = c("H", "Hex")
  to = "concrete"
  expect_snapshot(convert_mono_type(before, to), error = TRUE)
})


test_that("convert mono types with multiple `to`", {
  expect_error(
    convert_mono_type(c("Hex", "Hex"), to = c("simple", "simple")),
    "Only one `to` mono type can be specified"
  )
})


patrick::with_parameters_test_that("deciding glycan mono type works", {
    expect_equal(decide_glycan_mono_type(glycan), type)
  },
  glycan = list(
    n_glycan_core(mono_type = "concrete"),
    n_glycan_core(mono_type = "generic"),
    n_glycan_core(mono_type = "simple")
  ),
  type = c("concrete", "generic", "simple")
)


test_that("deciding mono types vectorized", {
  mono = c("Gal", "Hex", "H")
  expected = c("concrete", "generic", "simple")
  expect_equal(decide_mono_type(mono), expected)
})


test_that("deciding mono types fails for multiple monos", {
  mono = c("bad1", "bad2", "bad3")
  expect_snapshot(decide_mono_type(mono), error = TRUE)
})


test_that("ensure_glycan_mono_type works for the same mono type", {
  glycan <- o_glycan_core_1(mono_type = "concrete")
  glycan <- ensure_glycan_mono_type(glycan, "concrete")
  expect_equal(decide_glycan_mono_type(glycan), "concrete")
})


test_that("ensure_glycan_mono_type works for different mono types", {
  glycan <- o_glycan_core_1(mono_type = "concrete")
  glycan <- ensure_glycan_mono_type(glycan, "simple")
  expect_equal(decide_glycan_mono_type(glycan), "simple")
})


test_that("ensure_glycan_mono_type fails for bad conversion", {
  glycan <- o_glycan_core_1(mono_type = "simple")
  expect_error(ensure_glycan_mono_type(glycan, "concrete"))
})
