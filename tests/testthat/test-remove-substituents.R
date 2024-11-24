test_that("remove_substituents works for NE glycans", {
  glycan <- o_glycan_core_2(mode = "ne", linkage = TRUE)
  igraph::V(glycan)$sub[[1]] <- "3Me"
  glycan <- remove_substituents(glycan)
  expect_equal(igraph::V(glycan)$sub[[1]], "")
})


test_that("remove_substituents works for DN glycans", {
  glycan <- o_glycan_core_2(mode = "dn", linkage = TRUE)
  stopifnot(igraph::V(glycan)$type[[1]] == "mono")
  igraph::V(glycan)$sub[[1]] <- "3Me"
  glycan <- remove_substituents(glycan)
  expect_equal(igraph::V(glycan)$sub[[1]], "")
})
