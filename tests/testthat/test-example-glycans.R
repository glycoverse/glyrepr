test_that("N-glycan core NE graph", {
  glycan <- n_glycan_core(mode = "ne")
  expect_s3_class(glycan, c("ne_glycan_graph", "glycan_graph", "igraph"))
  expect_snapshot(str(igraph::edge_attr(glycan)))
  expect_snapshot(str(igraph::vertex_attr(glycan)))
})


test_that("N-glycan core DN graph", {
  glycan <- n_glycan_core(mode = "dn")
  expect_s3_class(glycan, c("dn_glycan_graph", "glycan_graph", "igraph"))
  expect_snapshot(str(igraph::edge_attr(glycan)))
  expect_snapshot(str(igraph::vertex_attr(glycan)))
})


test_that("N-glycan core graph without linkages", {
  ne_glycan <- n_glycan_core(mode = "ne", linkage = FALSE)
  dn_glycan <- n_glycan_core(mode = "dn", linkage = FALSE)

  expect_true(all(is.na(igraph::E(ne_glycan)$linkage)))
  expect_true(all(is.na(igraph::V(dn_glycan)$linkage)))
})
