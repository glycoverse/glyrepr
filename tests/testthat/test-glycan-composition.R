test_that("get the compositiona for a NE glycan graph", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::E(graph)$linkage <- NA_character_
  glycan <- new_ne_glycan_graph(graph)

  comp <- get_composition(glycan)

  expect_equal(sort(comp), sort(c(Glc = 2L, Gal = 1L)))
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


test_that("counting monos works for NE glycans", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::E(graph)$linkage <- NA_character_
  glycan <- new_ne_glycan_graph(graph)

  count <- count_monos(glycan)

  expect_equal(count, 3L)
})


test_that("counting monos works for DN glycans", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$type <- c("mono", "linkage", "mono")
  igraph::V(graph)$mono <- c("Glc", NA, "Glc")
  igraph::V(graph)$linkage <- c(NA, "b1-4", NA)
  glycan <- new_dn_glycan_graph(graph)

  count <- count_monos(glycan)

  expect_equal(count, 2L)
})
