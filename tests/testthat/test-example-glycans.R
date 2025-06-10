test_that("N-glycan core structure", {
  skip_on_old_win()
  glycan <- n_glycan_core()
  expect_s3_class(glycan, c("glyrepr_structure", "vctrs_rcrd", "vctrs_vctr"))
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("N-glycan core structure without linkages", {
  skip_on_old_win()
  glycan <- n_glycan_core(linkage = FALSE)
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("N-glycan core structure with generic monosaccharides", {
  skip_on_old_win()
  glycan <- n_glycan_core(mono_type = "generic")
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("O-glycan core 1", {
  skip_on_old_win()
  glycan <- o_glycan_core_1()
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("O-glycan core 2", {
  skip_on_old_win()
  glycan <- o_glycan_core_2()
  expect_snapshot(print(glycan, verbose = TRUE))
})
