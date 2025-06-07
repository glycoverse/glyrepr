test_that("N-glycan core graph", {
  skip_on_old_win()
  glycan <- n_glycan_core()
  expect_s3_class(glycan, c("glycan_graph", "igraph"))
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_glycan_graph(glycan))
})


test_that("N-glycan core graph without linkages", {
  skip_on_old_win()
  glycan <- n_glycan_core(linkage = FALSE)
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_glycan_graph(glycan))
})


test_that("N-glycan core graph with simple monosaccharides", {
  skip_on_old_win()
  glycan <- n_glycan_core(mono_type = "simple")
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_glycan_graph(glycan))
})



test_that("N-glycan core graph with generic monosaccharides", {
  skip_on_old_win()
  glycan <- n_glycan_core(mono_type = "generic")
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_glycan_graph(glycan))
})


test_that("O-glycan core 1", {
  skip_on_old_win()
  glycan <- o_glycan_core_1()
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_glycan_graph(glycan))
})


test_that("O-glycan core 2", {
  skip_on_old_win()
  glycan <- o_glycan_core_2()
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_glycan_graph(glycan))
})
