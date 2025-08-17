test_that("as_glycan_structure.character parses simple IUPAC-condensed strings", {
  # Single monosaccharide
  glycan1 <- as_glycan_structure("Man(?1-")
  expect_s3_class(glycan1, "glyrepr_structure")
  expect_equal(length(glycan1), 1)
  expect_equal(structure_to_iupac(glycan1), "Man(?1-")

  # Two monosaccharides
  glycan2 <- as_glycan_structure("Gal(b1-3)GalNAc(?1-")
  expect_s3_class(glycan2, "glyrepr_structure")
  expect_equal(length(glycan2), 1)
  expect_equal(structure_to_iupac(glycan2), "Gal(b1-3)GalNAc(?1-")

  # With explicit anomer
  glycan3 <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
  expect_s3_class(glycan3, "glyrepr_structure")
  expect_equal(structure_to_iupac(glycan3), "Gal(b1-3)GalNAc(a1-")
})

test_that("as_glycan_structure.character parses branched structures", {
  # Simple branched structure
  iupac <- "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"
  glycan <- as_glycan_structure(iupac)
  expect_s3_class(glycan, "glyrepr_structure")
  expect_equal(structure_to_iupac(glycan), "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
})

test_that("as_glycan_structure.character handles substituents", {
  # With substituents
  iupac <- "Man3S(a1-2)Gal6Ac(?1-"
  glycan <- as_glycan_structure(iupac)
  expect_s3_class(glycan, "glyrepr_structure")

  # Check that substituents are preserved
  graph <- get_structure_graphs(glycan, return_list = FALSE)
  expect_equal(igraph::V(graph)$sub, c("6Ac", "3S"))
})

test_that("as_glycan_structure.character handles Neu5Ac correctly", {
  # Neu5Ac with explicit anomer
  glycan1 <- as_glycan_structure("Neu5Ac(?2-")
  expect_equal(structure_to_iupac(glycan1), "Neu5Ac(?2-")

  # Neu5Ac with explicit anomer
  glycan2 <- as_glycan_structure("Neu5Ac(a2-")
  expect_equal(structure_to_iupac(glycan2), "Neu5Ac(a2-")

  # Neu5Ac with substituent
  glycan3 <- as_glycan_structure("Neu5Ac9Ac(?2-")
  graph <- get_structure_graphs(glycan3, return_list = FALSE)
  expect_equal(igraph::V(graph)$mono, "Neu5Ac")
  expect_equal(igraph::V(graph)$sub, "9Ac")
})



test_that("as_glycan_structure.character handles unknown linkages", {
  # Unknown linkages
  iupac <- "Man(a1-?)Man(?1-3)Man(?1-"
  glycan <- as_glycan_structure(iupac)
  expect_s3_class(glycan, "glyrepr_structure")
  expect_equal(structure_to_iupac(glycan), "Man(a1-?)Man(?1-3)Man(?1-")
})

test_that("as_glycan_structure.character handles multiple linkages", {
  # Multiple linkages
  iupac <- "Neu5Ac(a2-3/6)Gal(?1-"
  glycan <- as_glycan_structure(iupac)
  expect_s3_class(glycan, "glyrepr_structure")
  expect_equal(structure_to_iupac(glycan), "Neu5Ac(a2-3/6)Gal(?1-")
})

test_that("as_glycan_structure.character works with vectors", {
  # Multiple IUPAC strings
  iupacs <- c("Man(?1-", "Gal(b1-3)GalNAc(?1-", "Neu5Ac(a2-")
  glycans <- as_glycan_structure(iupacs)
  expect_s3_class(glycans, "glyrepr_structure")
  expect_equal(length(glycans), 3)

  # Check each one
  expect_equal(structure_to_iupac(glycans)[1], "Man(?1-")
  expect_equal(structure_to_iupac(glycans)[2], "Gal(b1-3)GalNAc(?1-")
  expect_equal(structure_to_iupac(glycans)[3], "Neu5Ac(a2-")
})

test_that("as_glycan_structure.character handles complex O-glycan", {
  # Complex O-glycan structure
  iupac <- "Neu5Ac(a2-3)Gal(b1-4)[Fuc(a1-3)]GlcNAc(b1-6)[Neu5Ac(a2-3)Gal(b1-3)]GalNAc(?1-"
  glycan <- as_glycan_structure(iupac)
  expect_s3_class(glycan, "glyrepr_structure")
  expect_equal(length(glycan), 1)

  # Check the structure has the correct number of nodes
  graph <- get_structure_graphs(glycan, return_list = FALSE)
  expect_equal(igraph::vcount(graph), 7)  # 7 monosaccharides
})

test_that("as_glycan_structure.character error handling", {
  # Empty string
  expect_error(as_glycan_structure(""), "Cannot parse empty")

  # NA
  expect_error(as_glycan_structure(NA_character_), "Cannot parse empty")

  # Invalid format - unknown monosaccharide
  expect_error(as_glycan_structure("invalid_format"), "Could not parse")

  # Missing anomer information
  expect_error(as_glycan_structure("Man"), "Can't extract anomer information")
  expect_error(as_glycan_structure("Neu5Ac"), "Can't extract anomer information")
  expect_error(as_glycan_structure("Gal(b1-3)GalNAc"), "Can't extract anomer information")
})

test_that("as_glycan_structure.character round-trip consistency", {
  # Test round-trip: structure -> IUPAC -> structure
  original_structures <- c(
    o_glycan_core_1(),
    n_glycan_core()
  )
  
  for (i in seq_along(original_structures)) {
    # Convert to IUPAC
    iupac <- structure_to_iupac(glycan_structure(original_structures[[i]]))
    
    # Parse back
    parsed <- as_glycan_structure(iupac)
    
    # Should be identical
    expect_equal(structure_to_iupac(parsed), iupac)
  }
})

# Edge cases and boundary conditions ----------------------------------------

test_that("as_glycan_structure.character handles whitespace", {
  # Leading and trailing whitespace should be handled gracefully
  expect_error(as_glycan_structure(" Man "), "Could not parse")
  expect_error(as_glycan_structure("  "), "Cannot parse empty")
  expect_error(as_glycan_structure("\t\n"), "Cannot parse empty")
})

test_that("as_glycan_structure.character handles malformed brackets", {
  # Unmatched brackets
  expect_error(as_glycan_structure("Man[Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man]Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man[[Gal]]"), "Could not parse")
  expect_error(as_glycan_structure("Man(a1-3)[Gal"), "Could not parse")
})

test_that("as_glycan_structure.character handles malformed linkages", {
  # Empty parentheses
  expect_error(as_glycan_structure("Man()Gal"), "Could not parse")
  
  # Incomplete linkages
  expect_error(as_glycan_structure("Man(a1)Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man(-3)Gal"), "Could not parse")
  
  # Invalid characters in linkages
  expect_error(as_glycan_structure("Man(c1-3)Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man(a0-3)Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man(a3-3)Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man(a1-0)Gal"), "Could not parse")
  expect_error(as_glycan_structure("Man(a#-3)Gal"), "Could not parse")
})

test_that("as_glycan_structure.character handles invalid monosaccharide names", {
  # Too short
  expect_error(as_glycan_structure("M"), "Could not parse")
  expect_error(as_glycan_structure("G"), "Could not parse")
  
  # Numbers at start
  expect_error(as_glycan_structure("1Man"), "Could not parse")
  expect_error(as_glycan_structure("2Gal"), "Could not parse")
  
  # Special characters
  expect_error(as_glycan_structure("Man@"), "Could not parse")
  expect_error(as_glycan_structure("Gal#"), "Could not parse")
  expect_error(as_glycan_structure("Man-Gal"), "Could not parse")
  
  # Case sensitivity (monosaccharide names are case sensitive)
  expect_error(as_glycan_structure("man"), "Could not parse")
  expect_error(as_glycan_structure("MAN"), "Could not parse")
  expect_error(as_glycan_structure("gal"), "Could not parse")
})

test_that("as_glycan_structure.character handles complex Neu variants", {
  # Test special Neu variants from glyparse
  glycan1 <- as_glycan_structure("Neu4Ac5Ac(?2-")
  graph1 <- get_structure_graphs(glycan1, return_list = FALSE)
  expect_equal(igraph::V(graph1)$mono, "Neu5Ac")
  expect_equal(igraph::V(graph1)$sub, "4Ac")

  glycan2 <- as_glycan_structure("Neu4Ac5Gc(?2-")
  graph2 <- get_structure_graphs(glycan2, return_list = FALSE)
  expect_equal(igraph::V(graph2)$mono, "Neu5Gc")
  expect_equal(igraph::V(graph2)$sub, "4Ac")

  # Neu5Ac should not be split even though it matches substituent pattern
  glycan3 <- as_glycan_structure("Neu5Ac(?2-")
  graph3 <- get_structure_graphs(glycan3, return_list = FALSE)
  expect_equal(igraph::V(graph3)$mono, "Neu5Ac")
  expect_equal(igraph::V(graph3)$sub, "")
})

test_that("as_glycan_structure.character handles invalid substituents", {
  # Invalid substituent positions
  expect_error(as_glycan_structure("Man0Ac"), "Could not parse")
  expect_error(as_glycan_structure("ManAAc"), "Could not parse")
  expect_error(as_glycan_structure("Man1X"), "Could not parse")  # X is not a valid substituent
})

test_that("as_glycan_structure.character handles deeply nested structures", {
  # Very deep nesting
  deep_structure <- "Man(a1-2)[Man(a1-4)[Man(a1-4)[Man(a1-5)Man(a1-2)]Man(a1-3)]Man(a1-4)]Man(a1-"
  glycan <- as_glycan_structure(deep_structure)
  expect_s3_class(glycan, "glyrepr_structure")
  
  # Verify the graph has correct number of nodes
  graph <- get_structure_graphs(glycan, return_list = FALSE)
  expect_equal(igraph::vcount(graph), 8)
})

test_that("as_glycan_structure.character handles mixed valid/invalid in vectors", {
  # Vector with some valid and some invalid strings
  mixed_vector <- c("Man", "invalid_mono", "Gal")
  expect_error(as_glycan_structure(mixed_vector), "Could not parse")
  
  # All invalid
  all_invalid <- c("invalid1", "invalid2")
  expect_error(as_glycan_structure(all_invalid), "Could not parse")
})



test_that("as_glycan_structure.character handles extreme linkage cases", {
  # Maximum valid positions
  glycan1 <- as_glycan_structure("Man(a2-9)Gal(?1-")
  expect_s3_class(glycan1, "glyrepr_structure")

  glycan2 <- as_glycan_structure("Man(b1-1)Gal(?1-")
  expect_s3_class(glycan2, "glyrepr_structure")

  # Complex multi-position linkages
  glycan3 <- as_glycan_structure("Man(a1-2/3/4/5/6)Gal(?1-")
  expect_s3_class(glycan3, "glyrepr_structure")
  expect_equal(structure_to_iupac(glycan3), "Man(a1-2/3/4/5/6)Gal(?1-")
})

test_that("as_glycan_structure.character handles single character edge cases", {
  # Very short but invalid inputs
  expect_error(as_glycan_structure("("), "Could not parse")
  expect_error(as_glycan_structure(")"), "Could not parse")
  expect_error(as_glycan_structure("["), "Could not parse")
  expect_error(as_glycan_structure("]"), "Could not parse")
  expect_error(as_glycan_structure("-"), "Could not parse")
})

test_that("as_glycan_structure.character handles unusual anomer cases", {
  # Various anomer edge cases
  glycan1 <- as_glycan_structure("Neu5Ac(?2-")
  expect_equal(structure_to_iupac(glycan1), "Neu5Ac(?2-")
  
  glycan2 <- as_glycan_structure("Man(??-")
  expect_equal(structure_to_iupac(glycan2), "Man(??-")
  
  glycan3 <- as_glycan_structure("Fuc(?1-")
  expect_equal(structure_to_iupac(glycan3), "Fuc(?1-")
})

test_that("as_glycan_structure.character preserves complex substituent patterns", {
  # Test various substituent combinations
  glycan1 <- as_glycan_structure("Gal6S(b1-3)GlcNAc4S(?1-")
  graph1 <- get_structure_graphs(glycan1, return_list = FALSE)
  expect_equal(sort(igraph::V(graph1)$sub), sort(c("4S", "6S")))

  # Unknown position substituents
  glycan2 <- as_glycan_structure("Man?S(?1-")
  graph2 <- get_structure_graphs(glycan2, return_list = FALSE)
  expect_equal(igraph::V(graph2)$sub, "?S")
})

test_that("as_glycan_structure.character handles multiple substituents", {
  glycan <- as_glycan_structure("Glc3Me6S(a1-")
  graph <- get_structure_graphs(glycan, return_list = FALSE)
  expect_equal(igraph::V(graph)$sub, "3Me,6S")
})

test_that("Neu monosaccharides with 5Ac are correctly parsed as Neu5Ac", {
  # Test cases where 5Ac should result in Neu5Ac base monosaccharide
  expect_equal(.extract_substituent("Neu3Me5Ac"), c(mono = "Neu5Ac", sub = "3Me"))
  expect_equal(.extract_substituent("Neu4Ac5Ac"), c(mono = "Neu5Ac", sub = "4Ac"))
  expect_equal(.extract_substituent("Neu4Ac5Ac9Ac"), c(mono = "Neu5Ac", sub = "4Ac,9Ac"))
  expect_equal(.extract_substituent("Neu7S5Ac"), c(mono = "Neu5Ac", sub = "7S"))
})

test_that("Neu monosaccharides with 5Gc are correctly parsed as Neu5Gc", {
  # Test cases where 5Gc should result in Neu5Gc base monosaccharide
  expect_equal(.extract_substituent("Neu3Me5Gc"), c(mono = "Neu5Gc", sub = "3Me"))
  expect_equal(.extract_substituent("Neu4Ac5Gc"), c(mono = "Neu5Gc", sub = "4Ac"))
  expect_equal(.extract_substituent("Neu7S5Gc"), c(mono = "Neu5Gc", sub = "7S"))
})

test_that("Neu monosaccharides without 5Ac or 5Gc remain as Neu", {
  # Test cases where no 5Ac or 5Gc should result in Neu base monosaccharide
  expect_equal(.extract_substituent("Neu"), c(mono = "Neu", sub = ""))
  expect_equal(.extract_substituent("Neu7Ac"), c(mono = "Neu", sub = "7Ac"))
  expect_equal(.extract_substituent("Neu3Me7Ac"), c(mono = "Neu", sub = "3Me,7Ac"))
})

test_that("Neu5Ac and Neu5Gc exact matches work correctly", {
  # Test exact matches
  expect_equal(.extract_substituent("Neu5Ac"), c(mono = "Neu5Ac", sub = ""))
  expect_equal(.extract_substituent("Neu5Gc"), c(mono = "Neu5Gc", sub = ""))
})

test_that("Neu5Ac and Neu5Gc with additional substituents work correctly", {
  # Test Neu5Ac/Neu5Gc with additional substituents
  expect_equal(.extract_substituent("Neu5Ac9Ac"), c(mono = "Neu5Ac", sub = "9Ac"))
  expect_equal(.extract_substituent("Neu5Gc9Ac"), c(mono = "Neu5Gc", sub = "9Ac"))
  expect_equal(.extract_substituent("Neu5Ac7S9Ac"), c(mono = "Neu5Ac", sub = "7S,9Ac"))
})

test_that("Error is thrown for monosaccharides with both 5Ac and 5Gc", {
  # This should be an error case
  expect_error(.extract_substituent("Neu5Ac5Gc"), "cannot have both 5Ac and 5Gc")
  expect_error(.extract_substituent("Neu5Gc5Ac"), "cannot have both 5Ac and 5Gc")
})

test_that("Full IUPAC parsing works with corrected Neu substituents", {
  # Test full IUPAC parsing
  result1 <- suppressMessages(as_glycan_structure("Neu3Me5Ac(a2-3)Gal(b1-4)Glc(a1-"))
  expect_true(is_glycan_structure(result1))
  
  result2 <- suppressMessages(as_glycan_structure("Neu4Ac5Gc(a2-3)Gal(b1-4)Glc(a1-"))
  expect_true(is_glycan_structure(result2))
  
  # Check that the output contains the correct monosaccharide names
  expect_match(as.character(result1), "Neu5Ac3Me")
  expect_match(as.character(result2), "Neu5Gc4Ac")
})
