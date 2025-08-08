test_that("get_anomer works", {
  x <- n_glycan_core()
  expect_equal(get_anomer(x), "b1")
})
