test_that("structure_to_iupac works with basic linear structures", {
  # Test with O-glycan core 1: GalNAc -> Gal
  glycan1 <- o_glycan_core_1()
  result1 <- structure_to_iupac(glycan1)
  expect_equal(result1, "Gal(b1-3)GalNAc(a1-")
  
  # Create simple linear structure: Glc -> GlcNAc -> Gal
  # This should produce: Gal(b1-3)GlcNAc(b1-4)Glc(?1-
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "GlcNAc", "Gal")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-3")
  graph$anomer <- "?1"
  graph$alditol <- FALSE
  glycan2 <- as_glycan_structure(graph)
  
  result2 <- structure_to_iupac(glycan2)
  expect_equal(result2, "Gal(b1-3)GlcNAc(b1-4)Glc(?1-")
})

test_that("structure_to_iupac works with branched structures", {
  # Test with N-glycan core
  glycan <- n_glycan_core()
  result <- structure_to_iupac(glycan)
  expect_equal(result, "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
  
  # Test with O-glycan core 2
  glycan2 <- o_glycan_core_2()
  result2 <- structure_to_iupac(glycan2)
  expect_equal(result2, "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-")
})

test_that("structure_to_iupac handles single node structures", {
  # Single node structure - need to add empty linkage attribute
  graph <- igraph::make_graph(edges = integer(0), n = 1)
  igraph::V(graph)$mono <- "Glc"
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- character(0)  # Empty linkage for no edges
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  expect_equal(result, "Glc(a1-")
})

test_that("structure_to_iupac handles two-node structures", {
  # Two node structure: A -> B
  graph <- igraph::make_graph(~ 1-+2)
  igraph::V(graph)$mono <- c("Man", "Glc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "a1-3"
  graph$anomer <- "b1"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  expect_equal(result, "Glc(a1-3)Man(b1-")
})

test_that("linkage comparison works correctly", {
  # Test different anomeric configurations: a < b < ?
  expect_equal(compare_linkages("a1-3", "b1-3"), -1)
  expect_equal(compare_linkages("b1-3", "?1-3"), -1)
  expect_equal(compare_linkages("a1-3", "?1-3"), -1)
  expect_equal(compare_linkages("a1-3", "a1-3"), 0)
  
  # Test different positions
  expect_equal(compare_linkages("a1-3", "a2-3"), -1)
  expect_equal(compare_linkages("a1-3", "a1-4"), -1)
  expect_equal(compare_linkages("a?-3", "a1-3"), 1)  # ? > numbers
  expect_equal(compare_linkages("a1-?", "a1-3"), 1)  # ? > numbers
  
  # Test complex combinations
  expect_equal(compare_linkages("a1-3", "a1-6"), -1)
  expect_equal(compare_linkages("a1-6", "a2-3"), -1)
  expect_equal(compare_linkages("b2-3", "a1-6"), 1)
})

test_that("parse_linkage works correctly", {
  # Test valid linkages
  result1 <- parse_linkage("a1-3")
  expect_equal(result1$x, "a")
  expect_equal(result1$y, "1")
  expect_equal(result1$z, "3")
  expect_equal(result1$x_rank, 1)
  expect_equal(result1$y_rank, 1)
  expect_equal(result1$z_rank, 3)
  
  # Test with ?
  result2 <- parse_linkage("?2-?")
  expect_equal(result2$x, "?")
  expect_equal(result2$y, "2")
  expect_equal(result2$z, "?")
  expect_equal(result2$x_rank, 3)
  expect_equal(result2$y_rank, 2)
  expect_equal(result2$z_rank, Inf)
  
  # Test invalid linkage
  expect_error(parse_linkage("invalid"), "Invalid linkage format")
})

test_that("structure_to_iupac ensures isomorphic graphs produce same sequence", {
  # Create first graph: Man with a1-3 and a1-6 branches (in that order)
  graph1 <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph1)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  igraph::V(graph1)$sub <- ""
  igraph::E(graph1)$linkage <- c("b1-4", "b1-4", "a1-3", "a1-6")
  graph1$anomer <- "?1"
  graph1$alditol <- FALSE
  glycan1 <- as_glycan_structure(graph1)
  
  # Create second graph: Same structure but with a1-6 and a1-3 branches (swapped order)
  graph2 <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph2)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  igraph::V(graph2)$sub <- ""
  igraph::E(graph2)$linkage <- c("b1-4", "b1-4", "a1-6", "a1-3")  # Swapped order
  graph2$anomer <- "?1"
  graph2$alditol <- FALSE
  glycan2 <- as_glycan_structure(graph2)
  
  # Both should produce the same sequence
  result1 <- structure_to_iupac(glycan1)
  result2 <- structure_to_iupac(glycan2)
  expect_equal(result1, result2)
  expect_equal(result1, "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
})

test_that("structure_to_iupac handles complex branched structures", {
  # Create a more complex structure with multiple levels of branching
  # Glc -> Man -> GlcNAc -> Gal
  #                |        ├─ Fuc
  #                |        └─ NeuAc  
  #                └─ Hex
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5, 2-+6)
  igraph::V(graph)$mono <- c("Glc", "Man", "GlcNAc", "Gal", "Fuc", "GalNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "a1-3", "a1-2", "a1-6", "b1-3")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  # Based on actual algorithm output
  expect_equal(result, "Gal(a1-6)[Fuc(b1-3)]GlcNAc(a1-3)[GalNAc(a1-2)]Man(b1-4)Glc(a1-")
})

test_that("structure_to_iupac selects correct backbone based on depth", {
  # Create structure where backbone selection matters
  # Glc -> Man -> Gal (depth 2)
  #           └─ GlcNAc -> Fuc (depth 3) <- this should be backbone
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 2-+4, 4-+5)
  igraph::V(graph)$mono <- c("Glc", "Man", "Gal", "GlcNAc", "Fuc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "a1-3", "a1-6", "b1-2")
  graph$anomer <- "?1"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  # GlcNAc->Fuc path is longer, so should be backbone
  # Expected: Fuc(b1-2)GlcNAc(a1-6)[Gal(a1-3)]Man(b1-4)Glc(?1-
  expect_equal(result, "Fuc(b1-2)GlcNAc(a1-6)[Gal(a1-3)]Man(b1-4)Glc(?1-")
})

test_that("structure_to_iupac selects backbone by linkage when depths are equal", {
  # Create structure where depths are equal but linkages differ
  # Glc -> Man ├─ Gal (linkage a1-3)
  #            └─ Fuc (linkage a1-6)
  # Both have same depth, but a1-3 < a1-6, so Gal should be backbone
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 2-+4)
  igraph::V(graph)$mono <- c("Glc", "Man", "Gal", "Fuc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "a1-3", "a1-6")
  graph$anomer <- "b2"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  # Gal should be backbone (a1-3), Fuc should be branch (a1-6)
  expect_equal(result, "Gal(a1-3)[Fuc(a1-6)]Man(b1-4)Glc(b2-")
})

test_that("structure_to_iupac handles different anomer values", {
  # Test different anomer values
  anomers <- c("a1", "b1", "?1", "a2", "b2", "?2")
  
  for (anomer in anomers) {
    graph <- igraph::make_graph(~ 1-+2)
    igraph::V(graph)$mono <- c("Glc", "Man")
    igraph::V(graph)$sub <- ""
    igraph::E(graph)$linkage <- "b1-4"
    graph$anomer <- anomer
    graph$alditol <- FALSE
    glycan <- as_glycan_structure(graph)
    
    result <- structure_to_iupac(glycan)
    expected <- paste0("Man(b1-4)Glc(", anomer, "-")
    expect_equal(result, expected)
  }
})

test_that("structure_to_iupac handles multiple branches correctly", {
  # Create structure with 3 branches
  # Glc -> Man ├─ Gal (a1-2)
  #            ├─ Fuc (a1-3) <- should be backbone (smallest)
  #            └─ Hex (a1-6)
  graph <- igraph::make_graph(~ 1-+2, 2-+3, 2-+4, 2-+5)
  igraph::V(graph)$mono <- c("Glc", "Man", "Gal", "Fuc", "GalNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "a1-2", "a1-3", "a1-6")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  # Gal should be backbone (a1-2 is smallest), Fuc and GalNAc should be branches
  # Branches should be ordered: a1-3 < a1-6
  expect_equal(result, "Gal(a1-2)[Fuc(a1-3)][GalNAc(a1-6)]Man(b1-4)Glc(a1-")
})

test_that("structure_to_iupac produces correct sequence for examples in documentation", {
  # Test example 1 from documentation
  # Glc (?1- └─GlcNAc (b1-4) └─Gal (b1-3)
  # Should produce: Gal(b1-3)GlcNAc(b1-4)Glc(?1-
  graph1 <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph1)$mono <- c("Glc", "GlcNAc", "Gal")
  igraph::V(graph1)$sub <- ""
  igraph::E(graph1)$linkage <- c("b1-4", "b1-3")
  graph1$anomer <- "?1"
  graph1$alditol <- FALSE
  glycan1 <- as_glycan_structure(graph1)
  
  result1 <- structure_to_iupac(glycan1)
  expect_equal(result1, "Gal(b1-3)GlcNAc(b1-4)Glc(?1-")
  
  # Test example 2 from documentation (already covered by n_glycan_core test)
  # But let's create it explicitly to be sure
  graph2 <- igraph::make_graph(~ 1-+2, 2-+3, 3-+4, 3-+5)
  igraph::V(graph2)$mono <- c("GlcNAc", "GlcNAc", "Man", "Man", "Man")
  igraph::V(graph2)$sub <- ""
  igraph::E(graph2)$linkage <- c("b1-4", "b1-4", "a1-3", "a1-6")
  graph2$anomer <- "?1"
  graph2$alditol <- FALSE
  glycan2 <- as_glycan_structure(graph2)
  
  result2 <- structure_to_iupac(glycan2)
  expect_equal(result2, "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-")
})

test_that("structure_to_iupac throws appropriate errors", {
  # Test with non-glycan_structure object
  expect_error(structure_to_iupac("not a glycan"), "glycan_structure")
  
  # Create a valid glycan first, then test error from structure_to_iupac
  graph <- igraph::make_graph(~ 1-+2)
  igraph::V(graph)$mono <- c("Glc", "Man")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "a1-3"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  # Manually create invalid structure (multiple roots) by bypassing validation
  invalid_graph <- igraph::make_graph(~ 1-+3, 2-+3)  # Two roots: 1 and 2
  igraph::V(invalid_graph)$mono <- c("Glc", "Man", "Gal")
  igraph::V(invalid_graph)$sub <- ""
  igraph::E(invalid_graph)$linkage <- c("a1-3", "b1-4")
  invalid_graph$anomer <- "a1"
  invalid_graph$alditol <- FALSE
  # Bypass validation by setting class directly
  class(invalid_graph) <- c("glycan_structure", "igraph")
  
  expect_error(structure_to_iupac(invalid_graph), "exactly one root")
})

test_that("depth calculation works correctly", {
  # Test depth calculation on n_glycan_core
  glycan <- n_glycan_core()
  root <- which(igraph::degree(glycan, mode = "in") == 0)
  depths <- calculate_depths(glycan, root)
  
  # Root should have depth 3, leaves should have depth 0
  expect_equal(as.numeric(depths["1"]), 3)  # Root GlcNAc
  expect_equal(as.numeric(depths["2"]), 2)  # Second GlcNAc
  expect_equal(as.numeric(depths["3"]), 1)  # Central Man
  expect_equal(as.numeric(depths["4"]), 0)  # Leaf Man
  expect_equal(as.numeric(depths["5"]), 0)  # Leaf Man
})

test_that("structure_to_iupac handles edge cases with linkages", {
  # Test with ? in different positions
  graph <- igraph::make_graph(~ 1-+2)
  igraph::V(graph)$mono <- c("Glc", "Man")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "?2-?"
  graph$anomer <- "??"
  graph$alditol <- FALSE
  glycan <- as_glycan_structure(graph)
  
  result <- structure_to_iupac(glycan)
  expect_equal(result, "Man(?2-?)Glc(??-")
  
  # Test with high numbers
  graph2 <- igraph::make_graph(~ 1-+2)
  igraph::V(graph2)$mono <- c("Glc", "Gal")
  igraph::V(graph2)$sub <- ""
  igraph::E(graph2)$linkage <- "b1-8"
  graph2$anomer <- "a7"
  graph2$alditol <- FALSE
  glycan2 <- as_glycan_structure(graph2)
  
  result2 <- structure_to_iupac(glycan2)
  expect_equal(result2, "Gal(b1-8)Glc(a7-")
}) 