test_that("remove_substituents works for glycans", {
  glycan <- o_glycan_core_2(linkage = TRUE)
  igraph::V(glycan)$sub[[1]] <- "3Me"
  glycan <- remove_substituents(glycan)
  expect_equal(igraph::V(glycan)$sub[[1]], "")
})
