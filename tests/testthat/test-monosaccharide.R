patrick::with_parameters_test_that(
  "known monosaccharides",
  {
    expect_true(is_known_monosaccharide(mono))
  },
  mono = c("Gal", "Hex")
)


patrick::with_parameters_test_that(
  "unknown monosaccharides",
  {
    expect_false(is_known_monosaccharide(mono))
  },
  mono = c("X", "Hx", "Nac")
)


test_that("get all generic monosaccharides", {
  res <- available_monosaccharides("generic")
  expect_contains(res, c("Hex", "HexNAc", "dHex"))
  expect_true(length(unique(res)) == length(res))
})


test_that("get all concrete monosaccharides", {
  res <- available_monosaccharides("concrete")
  expect_contains(res, c("Gal", "Man", "GlcNAc"))
  expect_true(length(unique(res)) == length(res))
})


test_that("get all monosaccharides", {
  res <- available_monosaccharides()
  expect_contains(res, c("Hex", "Man"))
  expect_true(length(unique(res)) == length(res))
})


test_that("get_anomer_pos returns anomer positions for concrete monosaccharides", {
  monos <- c("Gal", "GalNAc", "Neu5Ac", "Kdo", "Fru")

  expect_equal(get_anomer_pos(monos), c(1L, 1L, 2L, 2L, 2L))
})


test_that("get_anomer_pos returns anomer positions for generic monosaccharides", {
  monos <- c("Hex", "HexNAc", "NeuAc", "gKdo")

  expect_equal(get_anomer_pos(monos), c(1L, 1L, 2L, 2L))
})


test_that("get_anomer_pos supports mixed generic and concrete monosaccharides", {
  monos <- c("Gal", "Hex", "Neu5Ac", "NeuAc")

  expect_equal(get_anomer_pos(monos), c(1L, 1L, 2L, 2L))
})


test_that("get_anomer_pos preserves names", {
  monos <- c(a = "Gal", b = "Neu5Ac")

  expect_equal(names(get_anomer_pos(monos)), names(monos))
})


test_that("get_anomer_pos rejects unknown monosaccharides", {
  expect_error(
    get_anomer_pos("X"),
    regexp = "known monosaccharide"
  )
})
