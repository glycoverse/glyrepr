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
  expect_error(glycan_composition(list(c(Hex = 1, HexNAc = 1))), "named integer")
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

# This test is replaced by tests for NA handling in character casting above

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

# Tests for NA handling in constructor ---------------------------------

test_that("glycan_composition accepts NULL to create NA", {
  comp <- glycan_composition(NULL)
  expect_equal(length(comp), 1)
  expect_true(is.na(comp))
})

test_that("glycan_composition accepts NA to create NA", {
  comp <- glycan_composition(NA)
  expect_equal(length(comp), 1)
  expect_true(is.na(comp))
})

test_that("glycan_composition handles mixed valid and NA", {
  comp <- glycan_composition(c(Hex = 5), NULL, c(Hex = 3))
  expect_equal(length(comp), 3)
  expect_false(is.na(comp[1]))
  expect_true(is.na(comp[2]))
  expect_false(is.na(comp[3]))
})

test_that("glycan_composition rejects named NA (likely typo)", {
  # c(Hex = NA) is likely a typo for c(Hex = 5), should error
  expect_error(glycan_composition(c(Hex = NA)), "positive")
})

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

# Tests for NA handling in vec_restore --------------------------------

test_that("vec_restore skips NA in type checking", {
  comp1 <- glycan_composition(c(Hex = 5, HexNAc = 2))
  comp2 <- glycan_composition(NULL)

  # This should not error - NA should be skipped in type check
  combined <- c(comp1, comp2)
  expect_equal(length(combined), 2)
  expect_false(is.na(combined[1]))
  expect_true(is.na(combined[2]))
})

test_that("combining with NA at beginning works", {
  comp1 <- glycan_composition(c(Hex = 5, HexNAc = 2))
  combined <- c(NA, comp1)
  expect_equal(length(combined), 2)
  expect_true(is.na(combined[1]))
  expect_false(is.na(combined[2]))
})

# Tests for NA handling in character casting --------------------------------

test_that("as_glycan_composition handles NA characters", {
  chars <- c("Hex(5)HexNAc(2)", NA)
  comp <- as_glycan_composition(chars)
  expect_equal(length(comp), 2)
  expect_false(is.na(comp[1]))
  expect_true(is.na(comp[2]))
})

test_that("as_glycan_composition handles NA at beginning", {
  chars <- c(NA, "Hex(5)HexNAc(2)")
  comp <- as_glycan_composition(chars)
  expect_equal(length(comp), 2)
  expect_true(is.na(comp[1]))
  expect_false(is.na(comp[2]))
})

test_that("as_glycan_composition handles all NA characters", {
  chars <- c(NA, NA)
  comp <- as_glycan_composition(chars)
  expect_equal(length(comp), 2)
  expect_true(is.na(comp[1]))
  expect_true(is.na(comp[2]))
})

test_that("as_glycan_composition reports correct index with NA and invalid chars", {
  # Regression test: error message should report original indices, not filtered indices
  # When mixing NA and invalid strings, index 3 (not 2) should be reported
  chars <- c("Hex(5)", NA, "invalid")
  expect_error(
    as_glycan_composition(chars),
    "at index 3"
  )
})

test_that("as_glycan_composition reports correct index with NA at beginning and invalid", {
  chars <- c(NA, NA, "HexNAc(1)", "also_invalid")
  expect_error(
    as_glycan_composition(chars),
    "at index 4"
  )
})

# Tests for NA handling in format -----------------------------------------

test_that("format shows <NA> for NA compositions", {
  comp <- c(glycan_composition(c(Hex = 5)), NA)
  expect_equal(format(comp), c("Hex(5)", "<NA>"))
})

test_that("as.character shows <NA> for NA compositions", {
  comp <- c(glycan_composition(c(Hex = 5)), NA)
  expect_equal(as.character(comp), c("Hex(5)", "<NA>"))
})

test_that("print handles NA compositions", {
  comp <- c(glycan_composition(c(Hex = 5)), NA, glycan_composition(c(Hex = 3)))
  output <- capture.output(print(comp))
  expect_true(any(grepl("<NA>", output)))
})

test_that("tibble printing handles NA compositions", {
  comp <- c(glycan_composition(c(Hex = 5)), NA)
  tibble <- tibble::tibble(comp = comp, a = 1:2)
  output <- capture.output(print(tibble))
  expect_true(any(grepl("<NA>", output)))
})

# ===== NA Support Comprehensive Tests =====

test_that("is.na returns correct logical for compositions with NA", {
  comp <- glycan_composition(c(Hex = 5))
  expect_equal(is.na(comp), FALSE)

  comp_na <- c(glycan_composition(c(Hex = 5)), NA)
  expect_equal(is.na(comp_na), c(FALSE, TRUE))

  comp_all_na <- glycan_composition(NA_integer_, NA_integer_)
  expect_equal(is.na(comp_all_na), c(TRUE, TRUE))
})

test_that("as.logical preserves NA semantics", {
  # NA compositions should become NA, not TRUE
  comp <- glycan_composition(c(Hex = 5))
  expect_equal(as.logical(comp), FALSE)

  comp_na <- c(glycan_composition(c(Hex = 5)), NA)
  expect_equal(as.logical(comp_na), c(FALSE, NA))

  comp_all_na <- glycan_composition(NA_integer_, NA_integer_)
  expect_equal(as.logical(comp_all_na), c(NA, NA))
})

test_that("anyNA detects NA compositions", {
  comp <- glycan_composition(c(Hex = 5))
  expect_false(anyNA(comp))

  comp_na <- c(glycan_composition(c(Hex = 5)), NA)
  expect_true(anyNA(comp_na))
})

test_that("rep handles compositions with NA", {
  comp <- c(glycan_composition(c(Hex = 5)), NA)
  repeated <- rep(comp, 2)
  expect_equal(length(repeated), 4)
  expect_equal(is.na(repeated), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("rep handles compositions with NA at beginning", {
  # Note: c(NA, comp) doesn't preserve type in base R, use vctrs::vec_c for reliable behavior
  comp <- glycan_composition(c(Hex = 5))
  comps <- vctrs::vec_c(NA, comp)
  repeated <- rep(comps, 2)
  expect_equal(length(repeated), 4)
  expect_equal(is.na(repeated), c(TRUE, FALSE, TRUE, FALSE))
})

test_that("subsetting preserves NA", {
  comp <- c(glycan_composition(c(Hex = 5)), NA, glycan_composition(c(Hex = 3)))

  expect_true(is.na(comp[2]))
  expect_equal(is.na(comp[c(1, 3)]), c(FALSE, FALSE))
  expect_equal(length(comp[c(1, 3)]), 2)
})

test_that("as_glycan_composition handles list with NULL", {
  result <- as_glycan_composition(list(c(Hex = 5), NULL))
  expect_equal(length(result), 2)
  expect_false(is.na(result[1]))
  expect_true(is.na(result[2]))
})

test_that("combining multiple compositions with NA", {
  comp1 <- glycan_composition(c(Hex = 5, HexNAc = 2))
  comp2 <- glycan_composition(c(Hex = 3, HexNAc = 1))

  # Note: c(NA, comp) doesn't preserve type, so we use composition first
  combined <- c(comp1, NA, comp2, NA)
  expect_equal(length(combined), 4)
  expect_false(is.na(combined[1]))
  expect_true(is.na(combined[2]))
  expect_false(is.na(combined[3]))
  expect_true(is.na(combined[4]))
})

test_that("empty composition vector with NA works", {
  comps <- c(NA_character_, NA)
  result <- as_glycan_composition(comps)
  expect_equal(length(result), 2)
  expect_true(is.na(result[1]))
  expect_true(is.na(result[2]))
})

test_that("format preserves order with mixed NA and valid", {
  comp <- glycan_composition(c(Hex = 2, HexNAc = 1))
  # Use composition first, then NA - c(NA, comp) doesn't preserve type in base R
  comps <- c(comp, NA, glycan_composition(c(Hex = 3)), NA)
  formatted <- format(comps)
  # Valid elements show as composition string, NA elements show as "<NA>"
  expect_equal(formatted[1], "Hex(2)HexNAc(1)")
  expect_equal(formatted[2], "<NA>")
  expect_equal(formatted[3], "Hex(3)")
  expect_equal(formatted[4], "<NA>")
})

test_that("NA compositions are correctly restored after subsetting", {
  comps <- c(glycan_composition(c(Hex = 5)), NA, glycan_composition(c(Hex = 3)))
  subset <- comps[c(1, 3)]
  expect_equal(length(subset), 2)
  expect_false(is.na(subset[1]))
  expect_false(is.na(subset[2]))
})

test_that("combining concrete compositions with NA preserves type", {
  comps <- c(glycan_composition(c(Glc = 5, Gal = 2)), NA)
  expect_s3_class(comps, "glyrepr_composition")
  expect_false(is.na(comps[1]))
  expect_true(is.na(comps[2]))
  expect_equal(format(comps[1]), "Glc(5)Gal(2)")
})
