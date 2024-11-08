patrick::with_parameters_test_that("known monosaccharides", {
  expect_true(is_known_monosaccharide(mono))
}, mono = c("Gal", "Hex", "H"))


patrick::with_parameters_test_that("unknown monosaccharides", {
  expect_false(is_known_monosaccharide(mono))
}, mono = c("X", "Hx", "Nac"))
