# Tests for count_mono function

# Test helper function to create simple glycan graph
create_test_glycan_graph <- function(mono_names, linkages = NULL, anomer = "?1") {
  n_nodes <- length(mono_names)
  if (n_nodes == 1) {
    graph <- igraph::make_empty_graph(n = 1)
    linkages <- character(0)
  } else {
    # Create linear chain: 1-+2-+3-+...
    edges <- c()
    for (i in 1:(n_nodes - 1)) {
      edges <- c(edges, i, i + 1)
    }
    graph <- igraph::make_graph(edges = edges, directed = TRUE)
    if (is.null(linkages)) {
      linkages <- rep("b1-4", n_nodes - 1)
    }
  }
  
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))
  igraph::V(graph)$mono <- mono_names
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- linkages
  graph$anomer <- anomer
  graph
}

# Tests for count_mono with glycan compositions ---------------------------

test_that("count_mono works with compositions", {
  # Test generic composition
  comp <- glycan_composition(c(Hex = 3, HexNAc = 2, dHex = 1))
  expect_equal(count_mono(comp, "Hex"), 3)
  expect_equal(count_mono(comp, "HexNAc"), 2)
  expect_equal(count_mono(comp, "HexA"), 0)  # Not present
  
  # Test concrete composition
  comp2 <- glycan_composition(c(Glc = 2, Gal = 1, GlcNAc = 2))
  expect_equal(count_mono(comp2, "Glc"), 2)
  expect_equal(count_mono(comp2, "Man"), 0)  # Not present
})

test_that("count_mono works for compositions unable to be converted to generic", {
  comp <- glycan_composition(c(GlcNAc = 2, Kdn = 1))
  expect_equal(count_mono(comp, "HexNAc"), 2)
  expect_equal(count_mono(comp, "Kdn"), 1)
})

test_that("count_mono works for compositions with special monosaccharides", {
  comp <- glycan_composition(c(gNeu = 1, Hex = 1))
  expect_equal(count_mono(comp, "gNeu"), 1)
  expect_equal(count_mono(comp, "Hex"), 1)
  expect_equal(count_mono(comp, "Glc"), NA_integer_)
})

test_that("count_mono works when counting generic in concrete compositions", {
  # When mono is generic, it should count all matching concrete monos
  comp <- glycan_composition(c(Glc = 2, Gal = 1, Man = 1, GlcNAc = 2, GalNAc = 1))
  
  # Hex should count Glc, Gal, Man (all hexoses)
  expect_equal(count_mono(comp, "Hex"), 4)  # 2 Glc + 1 Gal + 1 Man
  
  # HexNAc should count GlcNAc, GalNAc
  expect_equal(count_mono(comp, "HexNAc"), 3)  # 2 GlcNAc + 1 GalNAc
})

test_that("count_mono returns NA when counting concrete in generic compositions", {
  comp <- glycan_composition(c(Hex = 2, HexNAc = 1))
  expect_equal(count_mono(comp, "GalNAc"), NA_integer_)
})

test_that("count_mono works with multiple compositions", {
  # Test with multiple compositions in a vector
  comp_vec <- glycan_composition(
    c(Hex = 5, HexNAc = 2),
    c(Hex = 3, HexNAc = 1, dHex = 1),
    c(HexNAc = 4)
  )
  
  expect_equal(count_mono(comp_vec, "Hex"), c(5, 3, 0))
  expect_equal(count_mono(comp_vec, "HexNAc"), c(2, 1, 4))
  expect_equal(count_mono(comp_vec, "dHex"), c(0, 1, 0))
})

test_that("count_mono works with `mono` as NULL", {
  comp <- glycan_composition(c(Man = 5, GlcNAc = 2), c(Gal = 1, Man = 1, GalNAc = 1))
  expect_equal(count_mono(comp), c(7L, 3L))
})

test_that("`include_subs` works when `mono` is NULL", {
  comp <- glycan_composition(c(Glc = 1, S = 1))
  expect_equal(count_mono(comp), 1L)
  expect_equal(count_mono(comp, include_subs = TRUE), 2L)
})

test_that("count_mono works for substituents", {
  comp <- glycan_composition(c(Glc = 1, S = 1), c(Glc = 1))
  expect_equal(count_mono(comp, "S"), c(1L, 0L))
})

# Tests for count_mono with glycan structures ----------------------------

test_that("count_mono works with glycan structures", {
  # Test with simple structure
  graph <- create_test_glycan_graph(c("GlcNAc", "Gal", "Glc"))
  struct <- glycan_structure(graph)
  expect_equal(count_mono(struct, "GlcNAc"), 1)
  expect_equal(count_mono(struct, "Gal"), 1)
  expect_equal(count_mono(struct, "Man"), 0)  # Not present
  
  # Test with generic monos in structures  
  graph2 <- create_test_glycan_graph(c("Glc", "Gal", "Man", "GlcNAc"))
  struct2 <- glycan_structure(graph2)
  expect_equal(count_mono(struct2, "Hex"), 3)  # Count Glc, Gal, Man
  expect_equal(count_mono(struct2, "HexNAc"), 1)  # Count GlcNAc
})

test_that("count_mono works with multiple structures", {
  # Test with N-glycan and O-glycan cores
  n_glycan <- n_glycan_core()
  o_glycan <- o_glycan_core_1()
  struct_vec <- c(n_glycan, o_glycan)
  
  # N-glycan core has: 2 GlcNAc, 3 Man; O-glycan core has: 1 GalNAc, 1 Gal
  expect_equal(count_mono(struct_vec, "GlcNAc"), c(2, 0))
  expect_equal(count_mono(struct_vec, "GalNAc"), c(0, 1))
})

# Tests for parameter validation and edge cases ---------------------------

test_that("count_mono validates parameters and handles edge cases", {
  comp <- glycan_composition(c(Hex = 2, HexNAc = 1))
  
  # Parameter validation
  expect_error(count_mono(comp, "Unknown"), "must be a known monosaccharide")
  expect_error(count_mono(comp, c("Hex", "HexNAc")), "Must have length 1")
  expect_error(count_mono(comp, 123), "Must be of type 'string'")
  
  # Edge cases
  empty_comp <- glycan_composition()
  expect_equal(length(count_mono(empty_comp, "Hex")), 0)
  
  # Return type and vector length
  result <- count_mono(comp, "Hex")
  expect_type(result, "integer")
  expect_equal(result, 2L)
}) 