test_that("empty composition is valid", {
  comp <- glycan_composition()
  expect_s3_class(comp, "glyrepr_composition")
})

test_that("generic composition is valid", {
  comp <- glycan_composition(c(Hex = 2, HexNAc = 1))
  expect_s3_class(comp, "glyrepr_composition")
})

test_that("concrete composition is valid", {
  comp <- glycan_composition(c(Glc = 2, Gal = 1))
  expect_s3_class(comp, "glyrepr_composition")
})

test_that("compositions can contain substituents", {
  comp <- glycan_composition(c(Glc = 1, S = 1))
  expect_equal(as.character(comp), "Glc(1)S(1)")
})

test_that("compositions are sorted correctly", {
  # Hex should come before HexNAc based on the monosaccharides tibble
  comp <- glycan_composition(c(HexNAc = 1, Hex = 2))
  data <- vctrs::vec_data(comp)
  comp_data <- vctrs::field(data, "data")[[1]]
  expect_equal(names(comp_data), c("Hex", "HexNAc"))
})

test_that("substituents are located after monosaccharides", {
  comp <- glycan_composition(c(S = 1, Gal = 1, Ac = 1, Glc = 1))
  # 1. "Ac" and "S" are sorted after "Glc" and "Gal"
  # 2. Order of "Ac" and "S" is according to `available_substituents()`
  expect_equal(as.character(comp), "Glc(1)Gal(1)Ac(1)S(1)")
})

test_that("mixed types within one composition throws error", {
  expect_error(glycan_composition(c(Hex = 1, Glc = 1)), "Must have only one type of monosaccharide")
})

test_that("mixed types within one composition vector throws error", {
  expect_error(glycan_composition(c(Hex = 1, HexNAc = 1), c(Glc = 1, Gal = 1)), "Must have only one type of monosaccharide")
})

test_that("unknown monosaccharides throw error", {
  expect_error(glycan_composition(c(Glc = 1, unknown = 1)), "Must have only known monosaccharides")
})

test_that("glycan_composition rejects wrong types", {
  expect_error(glycan_composition(list(c(Hex = 1, HexNAc = 1))), "Must be one or more named integer vectors")
})

test_that("glycan_composition rejects empty compositions", {
  expect_error(glycan_composition(integer(0)), "at least one residue")
})

test_that("as_glycan_composition works with list of named vectors", {
  vec_list <- list(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1))
  comp <- as_glycan_composition(vec_list)
  expected <- glycan_composition(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1))
  expect_equal(comp, expected)
})

test_that("as_glycan_composition works with a single named vector", {
  comp <- as_glycan_composition(c(Hex = 5, HexNAc = 2))
  expected <- glycan_composition(c(Hex = 5, HexNAc = 2))
  expect_equal(comp, expected)
  expect_equal(comp, as_glycan_composition(list(c(Hex = 5, HexNAc = 2))))
})

test_that("as_glycan_composition returns existing composition unchanged", {
  original <- glycan_composition(c(Hex = 5, HexNAc = 2))
  result <- as_glycan_composition(original)
  expect_identical(result, original)
})

test_that("as_glycan_composition works for a glycan structure", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  glycan <- glycan_structure(graph)

  comp <- as_glycan_composition(glycan)

  expect_s3_class(comp, "glyrepr_composition")
  expected_comp <- glycan_composition(c(Glc = 2L, Gal = 1L))
  expect_equal(comp, expected_comp)
})

test_that("as_glycan_composition works for a glycan structure with substituents", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "3Me")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  glycan <- glycan_structure(graph)

  comp <- as_glycan_composition(glycan)

  expect_s3_class(comp, "glyrepr_composition")
  expected_comp <- glycan_composition(c(Glc = 2L, Gal = 1L, Me = 1L))
  expect_equal(comp, expected_comp)
})

test_that("as_glycan_composition works for a glycan structure with multiple substituents", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "3Me,6S")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  glycan <- glycan_structure(graph)

  comp <- as_glycan_composition(glycan)

  expect_s3_class(comp, "glyrepr_composition")
  expected_comp <- glycan_composition(c(Glc = 2L, Gal = 1L, Me = 1L, S = 1L))
  expect_equal(comp, expected_comp)
})

# Tests for formatting ----------------------------------------------

test_that("format works correctly", {
  # generic monosaccharides
  comp2 <- glycan_composition(c(Hex = 2, HexNAc = 1))
  expect_equal(format(comp2), "Hex(2)HexNAc(1)")
  
  # concrete monosaccharides
  comp3 <- glycan_composition(c(Glc = 2, Gal = 1))
  expect_equal(format(comp3), "Glc(2)Gal(1)")
})

test_that("is_glycan_composition works correctly", {
  comp <- glycan_composition(c(Hex = 2, HexNAc = 1))
  expect_true(is_glycan_composition(comp))
  expect_false(is_glycan_composition(c(Hex = 2, HexNAc = 1)))
  expect_false(is_glycan_composition("Hex(2)HexNAc(1)"))
})

# Tests for c() function (vec_ptype2 and vec_cast methods) -----------------------

test_that("c() combines composition vectors correctly", {
  # Test basic combination of composition vectors
  comp1 <- glycan_composition(c(Hex = 5, HexNAc = 2))
  comp2 <- glycan_composition(c(Hex = 3, HexNAc = 1))
  
  # This should work without error
  combined <- c(comp1, comp2)
  
  expect_s3_class(combined, "glyrepr_composition")
  expect_equal(length(combined), 2)
  
  # Check that both compositions are preserved
  formatted <- format(combined)
  expect_equal(formatted[1], "Hex(5)HexNAc(2)")
  expect_equal(formatted[2], "Hex(3)HexNAc(1)")
})

test_that("c() handles multiple composition vectors", {
  comp1 <- glycan_composition(c(Hex = 2, HexNAc = 1))
  comp2 <- glycan_composition(c(Hex = 3, HexNAc = 2), c(Hex = 1, HexNAc = 1))
  comp3 <- glycan_composition(c(Hex = 4, HexNAc = 3))
  
  combined <- c(comp1, comp2, comp3)
  
  expect_s3_class(combined, "glyrepr_composition")
  expect_equal(length(combined), 4)  # 1 + 2 + 1 = 4 total compositions
  
  formatted <- format(combined)
  expect_equal(formatted[1], "Hex(2)HexNAc(1)")
  expect_equal(formatted[2], "Hex(3)HexNAc(2)") 
  expect_equal(formatted[3], "Hex(1)HexNAc(1)")
  expect_equal(formatted[4], "Hex(4)HexNAc(3)")
})

test_that("c() works with empty composition vectors", {
  comp1 <- glycan_composition()  # Empty composition vector
  comp2 <- glycan_composition(c(Hex = 2, HexNAc = 1))
  
  combined1 <- c(comp1, comp2)
  combined2 <- c(comp2, comp1)
  
  expect_s3_class(combined1, "glyrepr_composition")
  expect_s3_class(combined2, "glyrepr_composition")
  expect_equal(length(combined1), 1)
  expect_equal(length(combined2), 1)
  expect_equal(format(combined1), "Hex(2)HexNAc(1)")
  expect_equal(format(combined2), "Hex(2)HexNAc(1)")
})

test_that("c() throws error for different monosaccharide types", {
  # Test that combining generic and concrete compositions throws an error
  generic_comp <- glycan_composition(c(Hex = 1, HexNAc = 1))
  concrete_comp <- glycan_composition(c(Glc = 1, Gal = 1))

  expect_error(c(generic_comp, concrete_comp), "Can't combine")
})

test_that("c() throws error for different monosaccharide types with substituents", {
  # Test that combining generic and concrete compositions throws an error
  generic_comp <- glycan_composition(c(Hex = 1, HexNAc = 1, Me = 1))
  concrete_comp <- glycan_composition(c(Glc = 1, Gal = 1, S = 1))

  expect_error(c(generic_comp, concrete_comp), "Can't combine")
})

test_that("c() maintains proper ordering within compositions", {
  # Test that monosaccharide ordering is preserved during combination
  comp1 <- glycan_composition(c(GalNAc = 1, Gal = 2))  # Out of order input
  comp2 <- glycan_composition(c(GalNAc = 1, Glc = 2))  # Out of order input

  combined <- c(comp1, comp2)

  formatted <- format(combined)
  expect_equal(formatted[1], "Gal(2)GalNAc(1)")  # Should be reordered
  expect_equal(formatted[2], "Glc(2)GalNAc(1)")  # Should be reordered
})

# Tests for vector casting functionality ----------------------------------------

test_that("composition vectors can be subset and maintain structure", {
  comp <- glycan_composition(c(Hex = 1, HexNAc = 1), c(Hex = 2, HexNAc = 2), c(Hex = 3, HexNAc = 3))
  
  # Test subsetting
  subset1 <- comp[1]
  subset2 <- comp[c(1, 3)]
  
  expect_s3_class(subset1, "glyrepr_composition")
  expect_s3_class(subset2, "glyrepr_composition")
  expect_equal(length(subset1), 1)
  expect_equal(length(subset2), 2)
  
  expect_equal(format(subset1), "Hex(1)HexNAc(1)")
  expect_equal(format(subset2), c("Hex(1)HexNAc(1)", "Hex(3)HexNAc(3)"))
})

test_that("composition vectors can be repeated", {
  comp <- glycan_composition(c(Hex = 2, HexNAc = 1))
  
  # Test rep() function which uses vctrs casting methods
  repeated <- rep(comp, 3)
  
  expect_s3_class(repeated, "glyrepr_composition")
  expect_equal(length(repeated), 3)
  expect_equal(format(repeated), rep("Hex(2)HexNAc(1)", 3))
})

# Tests for printing performance ----------------------------------------------

test_that("large composition printing performance is acceptable", {
  # Create a large composition vector
  comps <- rep(list(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1, dHex = 1)), 500)
  large_comp <- do.call(glycan_composition, comps)
  
  # This should complete quickly without freezing
  start_time <- Sys.time()
  output <- capture.output(print(large_comp))
  end_time <- Sys.time()
  
  # Should take less than 1 second
  expect_true(as.numeric(end_time - start_time) < 1)
  
  # Should only show 10 compositions plus summary line
  expect_length(output, 12)  # 10 compositions + header + summary
  expect_match(output[12], "\\.\\.\\. \\(990 more not shown\\)")
})

test_that("truncation works in tibble for compositions", {
  comp <- glycan_composition(c(Hex = 1, HexNAc = 1), c(Hex = 2, HexNAc = 1), c(Hex = 3, dHex = 1))
  tibble <- tibble::tibble(comp = comp, a = 1)
  expect_snapshot(print(tibble, width = 30))
})

# Tests for character conversion ----------------------------------------------
test_that("as_glycan_composition works for compositions", {
  chars <- "Hex(5)HexNAc(2)"
  expected <- glycan_composition(c(Hex = 5, HexNAc = 2))
  expect_equal(as_glycan_composition(chars), expected)
})

test_that("as_glycan_composition works for simple compositions", {
  chars <- c("H5N2", "H5N4S1F1", "H5N4A1G1")
  expected <- glycan_composition(
    c(Hex = 5, HexNAc = 2),
    c(Hex = 5, HexNAc = 4, NeuAc = 1, dHex = 1),
    c(Hex = 5, HexNAc = 4, NeuAc = 1, NeuGc = 1)
  )
  expect_equal(as_glycan_composition(chars), expected)
})

test_that("as_glycan_composition rejects NA", {
  chars <- c("Hex(5)HexNAc(2)", NA)
  expect_error(as_glycan_composition(chars), "Cannot parse NA as glycan composition.")
})

test_that("as_glycan_composition works for empty characters", {
  expect_equal(as_glycan_composition(character()), glycan_composition())
})

test_that("as_glycan_composition rejects empty strings", {
  chars <- c("", "Hex(5)HexNAc(2)")
  expect_error(as_glycan_composition(chars))
})

test_that("as_glycan_composition raises error for illegal characters", {
  chars <- c("invalid", "Hex(5)HexNAc(2)")
  err_msg <- "Characters cannot be parsed as glycan compositions at index 1"
  expect_error(as_glycan_composition(chars), err_msg)
})

test_that("as_glycan_composition reorder residues", {
  chars <- c("Hex(2)HexNAc(1)", "HexNAc(1)Hex(2)")
  comp <- as_glycan_composition(chars)
  expected <- c("Hex(2)HexNAc(1)", "Hex(2)HexNAc(1)")
  expect_equal(format(comp), expected)
})

test_that("as.character works for compositions", {
  chars <- c("Hex(2)HexNAc(1)", "Hex(5)HexNAc(2)")
  comp <- as_glycan_composition(chars)
  expect_equal(as.character(comp), chars)
})

# Tests for NA helper functions ----------------------------------------------

test_that(".is_na_composition_elem detects NA correctly", {
  expect_true(.is_na_composition_elem(NULL))
  expect_false(.is_na_composition_elem(c(Hex = 5)))
  expect_false(.is_na_composition_elem(integer(0)))
})

test_that(".valid_composition_element validates correctly", {
  # Valid input
  expect_equal(.valid_composition_element(c(Hex = 5, HexNAc = 2)), c(Hex = 5L, HexNAc = 2L))

  # Invalid input - unnamed
  expect_error(.valid_composition_element(c(5, 2)), "named integer")

  # Invalid input - empty
  expect_error(.valid_composition_element(integer(0)), "at least one residue")

  # Invalid input - unknown monosaccharide
  expect_error(.valid_composition_element(c(Unknown = 5)), "known monosaccharides")

  # Invalid input - non-positive
  expect_error(.valid_composition_element(c(Hex = 0)), "positive numbers")
})
