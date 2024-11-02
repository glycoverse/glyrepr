test_that("N-glycan core NE graph", {
  skip_on_old_win()
  glycan <- n_glycan_core(mode = "ne")
  expect_s3_class(glycan, c("ne_glycan_graph", "glycan_graph", "igraph"))
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(glycan))
})


test_that("N-glycan core DN graph", {
  skip_on_old_win()
  glycan <- n_glycan_core(mode = "dn")
  expect_s3_class(glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_dn_glycan_graph(glycan))
})


test_that("N-glycan core NE graph without linkages", {
  skip_on_old_win()
  glycan <- n_glycan_core(mode = "ne", linkage = FALSE)
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(glycan))
})


test_that("N-glycan core NE graph with simple monosaccharides", {
  skip_on_old_win()
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(glycan))
})



test_that("N-glycan core NE graph with generic monosaccharides", {
  skip_on_old_win()
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(glycan))
})


test_that("O-glycan core 1", {
  skip_on_old_win()
  glycan <- o_glycan_core_1()
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(glycan))
})


test_that("O-glycan core 2", {
  skip_on_old_win()
  glycan <- o_glycan_core_2()
  expect_snapshot(print(glycan, verbose = TRUE))
  expect_no_error(validate_ne_glycan_graph(glycan))
})
