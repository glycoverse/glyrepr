test_that("monosaccharides with substituents are colored correctly", {
  expect_equal(get_mono_color("Neu5Ac9Ac"), "#A54399")  # Purple for Neu5Ac
  expect_equal(get_mono_color("Gal6S"), "#FFD400")      # Yellow for Gal
  expect_equal(get_mono_color("Man3Me"), "#00A651")     # Green for Man
  expect_equal(get_mono_color("Glc"), "#0072BC")        # Blue for Glc
})

test_that("special sialic acid cases are handled correctly", {
  # Test color assignment for special cases
  expect_equal(get_mono_color("Neu4Ac5Ac"), "#A54399")  # Purple for Neu5Ac
  expect_equal(get_mono_color("Neu4Ac5Gc"), "#8FCCE9")  # Light blue for Neu5Gc
})
