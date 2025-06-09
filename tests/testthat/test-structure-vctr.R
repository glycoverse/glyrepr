# Tests for structure_vctr functions

# Helper function to create a simple glycan structure
create_simple_glycan <- function(mono_names, linkages, anomer = "?1") {
  n_nodes <- length(mono_names)
  if (n_nodes == 1) {
    graph <- igraph::make_empty_graph(n = 1)
  } else {
    # Create linear chain: 1-+2-+3-+...
    edges <- c()
    for (i in 1:(n_nodes - 1)) {
      edges <- c(edges, i, i + 1)
    }
    graph <- igraph::make_graph(edges = edges, directed = TRUE)
  }
  
  igraph::V(graph)$mono <- mono_names
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- linkages
  graph$anomer <- anomer
  graph$alditol <- FALSE
  glycan_structure(graph)
}

# Tests for structure_vctr constructor --------------------------------------------

test_that("structure_vctr creates empty vector by default", {
  sv <- structure_vctr()
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(length(sv), 0)
  expect_length(attr(sv, "iupacs"), 0)
  expect_length(attr(sv, "structures"), 0)
})

test_that("structure_vctr works with single glycan structure", {
  glycan <- o_glycan_core_1()
  sv <- structure_vctr(list(glycan))
  
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(length(sv), 1)
  expect_length(attr(sv, "iupacs"), 1)
  expect_length(attr(sv, "structures"), 1)
  expect_equal(attr(sv, "iupacs")[[1]], "Gal(b1-3)GalNAc(a1-")
})

test_that("structure_vctr works with multiple different glycan structures", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- structure_vctr(list(glycan1, glycan2))
  
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "iupacs"), 2)
  expect_length(attr(sv, "structures"), 2)
})

test_that("structure_vctr removes duplicates based on IUPAC codes", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- o_glycan_core_1()  # Same structure
  sv <- structure_vctr(list(glycan1, glycan2))
  
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(length(sv), 2)  # Original vector has 2 elements
  expect_length(attr(sv, "iupacs"), 1)  # But only 1 unique structure
  expect_length(attr(sv, "structures"), 1)
})

test_that("structure_vctr handles structures with different IUPAC but same graph topology", {
  # Create two structures that are topologically the same but have different anomers
  glycan1 <- create_simple_glycan(c("Glc", "Gal"), "b1-4", "a1")
  glycan2 <- create_simple_glycan(c("Glc", "Gal"), "b1-4", "b1")
  
  sv <- structure_vctr(list(glycan1, glycan2))
  
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "iupacs"), 2)  # Different IUPAC codes
})

test_that("structure_vctr validates input", {
  expect_error(structure_vctr("not a list"), "Assertion on 'x' failed")
  expect_error(structure_vctr(list(1, 2)), "Must inherit from class 'glycan_structure'")
})

# Tests for new_structure_vctr -------------------------------------------------

test_that("new_structure_vctr creates valid structure_vctr object", {
  codes <- c("hash1", "hash2")
  iupacs <- c("seq1", "seq2")
  structures <- list(o_glycan_core_1(), n_glycan_core())
  names(structures) <- codes
  names(iupacs) <- codes
  
  sv <- new_structure_vctr(codes, iupacs, structures)
  
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(vctrs::vec_data(sv), codes)
  expect_equal(attr(sv, "iupacs"), iupacs)
  expect_equal(attr(sv, "structures"), structures)
})

test_that("new_structure_vctr works with empty inputs", {
  sv <- new_structure_vctr()
  
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(length(sv), 0)
})

# Tests for is_structure_vctr --------------------------------------------------

test_that("is_structure_vctr correctly identifies structure_vctr objects", {
  sv <- structure_vctr()
  expect_true(is_structure_vctr(sv))
  
  expect_false(is_structure_vctr(1))
  expect_false(is_structure_vctr("test"))
  expect_false(is_structure_vctr(list()))
  expect_false(is_structure_vctr(o_glycan_core_1()))
})

# Tests for format method -----------------------------------------------------

test_that("format.glyrepr_structure_vctr displays correct IUPAC sequences", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- structure_vctr(list(glycan1, glycan2))
  
  formatted <- format(sv)
  
  expect_type(formatted, "character")
  expect_length(formatted, 2)
  expect_true("Gal(b1-3)GalNAc(a1-" %in% formatted)
  expect_true("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-" %in% formatted)
})

test_that("format.glyrepr_structure_vctr handles empty vector", {
  sv <- structure_vctr()
  formatted <- format(sv)
  
  expect_type(formatted, "character")
  expect_length(formatted, 0)
})

test_that("format.glyrepr_structure_vctr handles duplicates correctly", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- o_glycan_core_1()
  sv <- structure_vctr(list(glycan1, glycan2))
  
  formatted <- format(sv)
  
  expect_length(formatted, 2)
  expect_equal(formatted[1], formatted[2])  # Both should show same IUPAC
  expect_equal(formatted[1], "Gal(b1-3)GalNAc(a1-")
})

# Tests for vctrs methods -----------------------------------------------------

test_that("vec_ptype_abbr.glyrepr_structure_vctr returns correct abbreviation", {
  sv <- structure_vctr()
  expect_equal(vctrs::vec_ptype_abbr(sv), "struc_vctr")
})

test_that("vec_ptype_full.glyrepr_structure_vctr returns correct full type", {
  sv <- structure_vctr()
  expect_equal(vctrs::vec_ptype_full(sv), "structure_vctr")
})

# Tests for print footer ------------------------------------------------------

test_that("obj_print_footer.glyrepr_structure_vctr displays unique count", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  glycan3 <- o_glycan_core_1()  # Duplicate
  sv <- structure_vctr(list(glycan1, glycan2, glycan3))
  
  output <- capture.output(obj_print_footer.glyrepr_structure_vctr(sv))
  
  expect_length(output, 1)
  expect_match(output, "# Unique structures: 2")
})

test_that("obj_print_footer.glyrepr_structure_vctr handles empty vector", {
  sv <- structure_vctr()
  
  output <- capture.output(obj_print_footer.glyrepr_structure_vctr(sv))
  
  expect_length(output, 1)
  expect_match(output, "# Unique structures: 0")
})

# Integration tests -----------------------------------------------------------

test_that("structure_vctr preserves glycan structures correctly", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- structure_vctr(list(glycan1, glycan2))
  
  # Check that we can retrieve the original structures
  structures <- attr(sv, "structures")
  expect_length(structures, 2)
  
  # The structures should be equivalent to originals (though order may differ)
  iupacs <- attr(sv, "iupacs")
  expected_iupacs <- c(structure_to_iupac(glycan1), structure_to_iupac(glycan2))
  expect_setequal(iupacs, expected_iupacs)
})

test_that("structure_vctr handles complex branched structures", {
  # Create a custom branched structure
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 2-+4)
  igraph::V(graph)$mono <- c("GlcNAc", "Man", "Gal", "Fuc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "a1-3", "a1-6")
  graph$anomer <- "b1"
  graph$alditol <- FALSE
  complex_glycan <- glycan_structure(graph)
  
  sv <- structure_vctr(list(complex_glycan))
  
  expect_s3_class(sv, "glyrepr_structure_vctr")
  expect_equal(length(sv), 1)
  
  # Check the IUPAC sequence is generated correctly
  formatted <- format(sv)
  expect_type(formatted, "character")
  expect_length(formatted, 1)
})

test_that("structure_vctr maintains hash uniqueness property", {
  # Create multiple identical structures
  glycans <- replicate(5, o_glycan_core_1(), simplify = FALSE)
  sv <- structure_vctr(glycans)
  
  expect_equal(length(sv), 5)  # Original length maintained
  expect_length(attr(sv, "structures"), 1)  # Only one unique structure stored
  expect_length(attr(sv, "iupacs"), 1)  # Only one unique IUPAC stored
  
  # All elements should format to the same string
  formatted <- format(sv)
  expect_true(all(formatted == formatted[1]))
}) 