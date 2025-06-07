test_that("get the composition for a glycan graph", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_glycan_graph(graph)

  comp <- get_composition(glycan)

  expect_s3_class(comp, "glyrepr_composition")
  expected_comp <- composition(c(Glc = 2L, Gal = 1L))
  expect_equal(comp, expected_comp)
})

test_that("empty composition is valid", {
  comp <- composition()
  expect_s3_class(comp, "glyrepr_composition")
})

test_that("format is correct", {
  # simple monosaccharides
  comp <- composition(c(H = 2, N = 1))
  expect_equal(format(comp), "H2N1")
  # generic monosaccharides
  comp <- composition(c(Hex = 2, HexNAc = 1))
  expect_equal(format(comp), "Hex(2)HexNAc(1)")
  # concrete monosaccharides
  comp <- composition(c(Glc = 2, Gal = 1))
  expect_equal(format(comp), "Glc(2)Gal(1)")
})

test_that("residues are reordered", {
  comp <- composition(c(N = 1, H = 2))
  expect_equal(format(comp), "H2N1")
})

test_that("mixed mono types are not allowed", {
  expect_error(composition(c(H = 1, Hex = 1)), "must have only one type of monosaccharide")
})

test_that("unknown monosaccharides are not allowed", {
  expect_error(composition(c(Glc = 1, unknown = 1)), "must have only known monosaccharides")
})
