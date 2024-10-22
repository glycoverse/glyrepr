test_that("N-glycan core NE graph", {
  glycan <- n_glycan_core(mode = "ne")
  expect_s3_class(glycan, c("ne_glycan_graph", "glycan_graph", "igraph"))
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("N-glycan core DN graph", {
  glycan <- n_glycan_core(mode = "dn")
  expect_s3_class(glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("N-glycan core NE graph without linkages", {
  glycan <- n_glycan_core(mode = "ne", linkage = FALSE)
  expect_snapshot(print(glycan, verbose = TRUE))
})


test_that("N-glycan core NE graph with simple monosaccharides", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "simple")
  expect_snapshot(print(glycan, verbose = TRUE))
})



test_that("N-glycan core NE graph with generic monosaccharides", {
  glycan <- n_glycan_core(mode = "ne", mono_type = "generic")
  expect_equal(igraph::V(glycan)$mono, c("HexNAc", "HexNAc", "Hex", "Hex", "Hex"))
})
