test_that("get the compositiona for a NE glycan graph", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::E(graph)$linkage <- NA_character_
  glycan <- new_ne_glycan_graph(graph)

  comp <- get_composition(glycan)

  expect_equal(sort(comp), sort(c(Glc = 2L, Gal = 1L)))
})


test_that("get compositions for a list of NE glycans", {
  graph1 <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph1)$mono <- c("Glc", "Gal", "Glc")
  igraph::E(graph1)$linkage <- NA_character_
  glycan1 <- new_ne_glycan_graph(graph1)

  graph2 <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph2)$mono <- c("Hex", "Hex", "HexNac")
  igraph::E(graph2)$linkage <- NA_character_
  glycan2 <- new_ne_glycan_graph(graph2)

  comps <- get_compositions(list(glycan1, glycan2))

  expect_equal(sort(comps[[1]]), sort(c(Glc = 2L, Gal = 1L)))
  expect_equal(sort(comps[[2]]), sort(c(Hex = 2L, HexNac = 1L)))
})


test_that("get the composition for a DN glycan graph", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Glc")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA)
  glycan <- new_dn_glycan_graph(graph)

  comp <- get_composition(glycan)

  expect_equal(comp, c(Glc = 2L))
})


test_that("get compositions for a list of DN glycans", {
  graph1 <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph1)$type <- c("mono", "linkage", "mono")
  igraph::V(graph1)$mono <- c("Glc", NA, "Glc")
  igraph::V(graph1)$linkage <- c(NA, "b1-4", NA)
  glycan1 <- new_dn_glycan_graph(graph1)

  graph2 <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph2)$type <- c("mono", "linkage", "mono")
  igraph::V(graph2)$mono <- c("Hex", NA, "HexNac")
  igraph::V(graph2)$linkage <- c(NA, "b1-4", NA)
  glycan2 <- new_dn_glycan_graph(graph2)

  comps <- get_compositions(list(glycan1, glycan2))

  expect_equal(comps[[1]], c(Glc = 2L))
  expect_equal(comps[[2]], c(Hex = 1L, HexNac = 1L))
})
