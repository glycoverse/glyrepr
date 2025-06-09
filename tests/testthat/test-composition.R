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

test_that("as_composition works for a glycan structure", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- glycan_structure(graph)

  comp <- as_composition(glycan)

  expect_s3_class(comp, "glyrepr_composition")
  expected_comp <- composition(c(Glc = 2L, Gal = 1L))
  expect_equal(comp, expected_comp)
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

# Tests for c() function (vec_ptype2 and vec_cast methods) -----------------------

test_that("c() combines composition vectors correctly", {
  # Test basic combination of composition vectors
  comp1 <- composition(c(H = 5, N = 2))
  comp2 <- composition(c(Hex = 3, HexNAc = 1))
  
  # This should work without error
  combined <- c(comp1, comp2)
  
  expect_s3_class(combined, "glyrepr_composition")
  expect_equal(length(combined), 2)
  
  # Check that both compositions are preserved
  formatted <- format(combined)
  expect_equal(formatted[1], "H5N2")
  expect_equal(formatted[2], "Hex(3)HexNAc(1)")
})

test_that("c() handles multiple composition vectors", {
  comp1 <- composition(c(H = 2, N = 1))
  comp2 <- composition(c(H = 3, N = 2), c(H = 1, N = 1))
  comp3 <- composition(c(H = 4, N = 3))
  
  combined <- c(comp1, comp2, comp3)
  
  expect_s3_class(combined, "glyrepr_composition")
  expect_equal(length(combined), 4)  # 1 + 2 + 1 = 4 total compositions
  
  formatted <- format(combined)
  expect_equal(formatted[1], "H2N1")
  expect_equal(formatted[2], "H3N2") 
  expect_equal(formatted[3], "H1N1")
  expect_equal(formatted[4], "H4N3")
})

test_that("c() works with empty composition vectors", {
  comp1 <- composition()  # Empty composition vector
  comp2 <- composition(c(H = 2, N = 1))
  
  combined1 <- c(comp1, comp2)
  combined2 <- c(comp2, comp1)
  
  expect_s3_class(combined1, "glyrepr_composition")
  expect_s3_class(combined2, "glyrepr_composition")
  expect_equal(length(combined1), 1)
  expect_equal(length(combined2), 1)
  expect_equal(format(combined1), "H2N1")
  expect_equal(format(combined2), "H2N1")
})

test_that("c() preserves different monosaccharide types", {
  # Test with different mono types - each should remain separate
  simple_comp <- composition(c(H = 2, N = 1))
  generic_comp <- composition(c(Hex = 1, HexNAc = 1))
  concrete_comp <- composition(c(Glc = 1, Gal = 1))
  
  combined <- c(simple_comp, generic_comp, concrete_comp)
  
  expect_s3_class(combined, "glyrepr_composition")
  expect_equal(length(combined), 3)
  
  formatted <- format(combined)
  expect_equal(formatted[1], "H2N1")
  expect_equal(formatted[2], "Hex(1)HexNAc(1)")
  expect_equal(formatted[3], "Glc(1)Gal(1)")  # Corrected order based on monosaccharide tibble
})

test_that("c() maintains proper ordering within compositions", {
  # Test that monosaccharide ordering is preserved during combination
  comp1 <- composition(c(N = 1, H = 2))  # Out of order input
  comp2 <- composition(c(HexNAc = 1, Hex = 2))  # Out of order input
  
  combined <- c(comp1, comp2)
  
  formatted <- format(combined)
  expect_equal(formatted[1], "H2N1")  # Should be reordered
  expect_equal(formatted[2], "Hex(2)HexNAc(1)")  # Should be reordered
})

# Tests for vector casting functionality ----------------------------------------

test_that("composition vectors can be subset and maintain structure", {
  comp <- composition(c(H = 1, N = 1), c(H = 2, N = 2), c(H = 3, N = 3))
  
  # Test subsetting
  subset1 <- comp[1]
  subset2 <- comp[c(1, 3)]
  
  expect_s3_class(subset1, "glyrepr_composition")
  expect_s3_class(subset2, "glyrepr_composition")
  expect_equal(length(subset1), 1)
  expect_equal(length(subset2), 2)
  
  expect_equal(format(subset1), "H1N1")
  expect_equal(format(subset2), c("H1N1", "H3N3"))
})

test_that("composition vectors can be repeated", {
  comp <- composition(c(H = 2, N = 1))
  
  # Test rep() function which uses vctrs casting methods
  repeated <- rep(comp, 3)
  
  expect_s3_class(repeated, "glyrepr_composition")
  expect_equal(length(repeated), 3)
  expect_equal(format(repeated), rep("H2N1", 3))
})
