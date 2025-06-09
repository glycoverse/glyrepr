test_that("remove_substituents works for glycans", {
  glycan_vec <- o_glycan_core_2(linkage = TRUE)
  glycan <- get_structure_graphs(glycan_vec, 1)  # Extract the igraph
  igraph::V(glycan)$sub[[1]] <- "3Me"
  
  # Create a new vectorized structure with the modified graph
  modified_glycan_vec <- glycan_structure(glycan)
  result <- remove_substituents(modified_glycan_vec)
  
  # Extract the result graph and check
  result_graph <- get_structure_graphs(result, 1)
  expect_equal(igraph::V(result_graph)$sub[[1]], "")
})
