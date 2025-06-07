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

test_that("as_composition works with named vectors", {
  # Test with simple monosaccharides
  vec <- c(H = 5, N = 2)
  comp <- as_composition(vec)
  expected <- composition(c(H = 5, N = 2))
  expect_equal(comp, expected)
  
  # Test with generic monosaccharides
  vec <- c(Hex = 2, HexNAc = 1)
  comp <- as_composition(vec)
  expected <- composition(c(Hex = 2, HexNAc = 1))
  expect_equal(comp, expected)
  
  # Test with concrete monosaccharides
  vec <- c(Glc = 2, Gal = 1)
  comp <- as_composition(vec)
  expected <- composition(c(Glc = 2, Gal = 1))
  expect_equal(comp, expected)
})

test_that("as_composition works with list of named vectors", {
  vec_list <- list(c(H = 5, N = 2), c(H = 3, N = 1))
  comp <- as_composition(vec_list)
  expected <- composition(c(H = 5, N = 2), c(H = 3, N = 1))
  expect_equal(comp, expected)
})

test_that("as_composition returns existing composition unchanged", {
  original <- composition(c(H = 5, N = 2))
  result <- as_composition(original)
  expect_identical(result, original)
})

test_that("as_composition handles edge cases", {
  # Empty named vector
  vec <- integer(0)
  names(vec) <- character(0)
  comp <- as_composition(vec)
  expected <- composition(vec)
  expect_equal(comp, expected)
  
  # Single element
  vec <- c(H = 1)
  comp <- as_composition(vec)
  expected <- composition(c(H = 1))
  expect_equal(comp, expected)
})

test_that("as_composition throws errors for invalid inputs", {
  # Unnamed vector
  expect_error(as_composition(c(1, 2, 3)), "Cannot convert")
  
  # Non-numeric vector
  expect_error(as_composition(c(a = "x", b = "y")), "Cannot convert")
  
  # Invalid object type
  expect_error(as_composition(data.frame(a = 1, b = 2)), "Cannot convert")
  
  # List with unnamed vectors
  expect_error(as_composition(list(c(1, 2))), "All elements in the list must be named vectors")
})

test_that("as_composition maintains sorting", {
  # Test that residues are properly sorted
  vec <- c(N = 1, H = 2)  # Out of order
  comp <- as_composition(vec)
  expect_equal(format(comp), "H2N1")  # Should be reordered
})
