patrick::with_parameters_test_that("known monosaccharides", {
  expect_true(is_known_monosaccharide(mono))
}, mono = c("Gal", "Hex", "H"))


patrick::with_parameters_test_that("unknown monosaccharides", {
  expect_false(is_known_monosaccharide(mono))
}, mono = c("X", "Hx", "Nac"))


test_that("get all simple monosaccharides", {
  expect_setequal(available_monosaccharides("simple"), c("H", "N", "F", "S", "A", "G"))
})


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
  expect_contains(res, c("H", "Hex", "Man"))
  expect_true(length(unique(res)) == length(res))
})
