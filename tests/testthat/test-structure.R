# Tests for glycan structure functions

good_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph
}


# Tests for glycan_structure --------------------------------------------------

test_that("glycan_structure works", {
  glycan <- glycan_structure(good_glycan_graph())
  expect_s3_class(glycan, c("glyrepr_structure"))
})


test_that("glycan_structure fails for invalid graphs", {
  bad_graph <- igraph::make_graph(~ 1-+2, 2-+3, 3-+1)
  expect_error(glycan_structure(bad_graph))
})


test_that("vertex names are added if missing", {
  graph <- good_glycan_graph()
  graph <- igraph::delete_vertex_attr(graph, "name")

  glycan_vec <- glycan_structure(graph)
  glycan <- get_structure_graphs(glycan_vec, 1)

  expect_true("name" %in% igraph::vertex_attr_names(glycan))
})


# Tests for glycan_structure validation during creation --------------------------

test_that("glycan_structure accepts valid graphs", {
  graph <- good_glycan_graph()
  expect_no_error(glycan_structure(graph))
})


test_that("glycan_structure rejects invalid graphs", {
  # Test undirected graph
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("Glc", "Gal")
  igraph::V(graph)$sub <- c("", "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "directed")
})


# Tests for ensure_name_vertex_attr ------------------------------------------

test_that("ensure_name_vertex_attr adds names when missing", {
  graph <- igraph::make_graph(~ 1-+2)
  graph <- igraph::delete_vertex_attr(graph, "name")
  result <- ensure_name_vertex_attr(graph)
  expect_true("name" %in% igraph::vertex_attr_names(result))
})


test_that("ensure_name_vertex_attr preserves existing names", {
  graph <- igraph::make_graph(~ A-+B)
  result <- ensure_name_vertex_attr(graph)
  expect_equal(igraph::V(result)$name, c("A", "B"))
})


test_that("glycan structure class", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  glycan <- glycan_structure(graph)
  expect_s3_class(glycan, c("glyrepr_structure"))
})


test_that("validating undirected graphs", {
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must be directed")
})


test_that("validating an in tree", {
  graph <- igraph::make_tree(3, children = 2, mode = "in")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must be an out tree")
})


test_that("validating graph without monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must have a vertex attribute 'mono'")
})


test_that("validating graph without substituent attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must have a vertex attribute 'sub'")
})


test_that("validating graph with NA in monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must have no NA in vertex attribute 'mono'")
})


test_that("validating graph with NA in substitude attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- c("", NA, "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must have no NA in vertex attribute 'sub'")
})


patrick::with_parameters_test_that("valid substituents", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- sub
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "b1"
  expect_no_error(glycan_structure(graph))
}, sub = c("6S", "9Ac", "2P", "?S"))


test_that("multiple substituents are supported", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "Glc"
  igraph::V(graph)$sub <- "3Me,4Ac"  # Multiple substituents
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "a1"
  expect_no_error(glycan_structure(graph))
})


test_that("multiple substituents must be sorted by position", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "Glc"
  igraph::V(graph)$sub <- "4Ac,3Me"  # Wrong order - should be "3Me,4Ac"
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Unknown substituent")
})


test_that("duplicate positions in substituents are not allowed", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "Glc"
  igraph::V(graph)$sub <- "3Me,3Ac"  # Same position (3) with different substituents
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Unknown substituent")
})


test_that("normalize_substituents works correctly", {
  expect_equal(normalize_substituents(""), "")
  expect_equal(normalize_substituents("6S"), "6S")
  expect_equal(normalize_substituents("4Ac,3Me"), "3Me,4Ac")
  expect_equal(normalize_substituents("6P,2S,4Ac"), "2S,4Ac,6P")
  expect_equal(normalize_substituents("?S,3Me"), "3Me,?S")
})


test_that("validating graph without linkage attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  expect_error(glycan_structure(graph), "Glycan structure must have an edge attribute 'linkage'")
})


test_that("validating one non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Fuc", "Bad")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"

  expect_error(glycan_structure(graph))
  err <- rlang::catch_cnd(glycan_structure(graph))
  # The error is now wrapped by purrr, so check the parent error
  expect_true(grepl("Unknown monosaccharide: Bad", err$parent$message))
  expect_equal(err$parent$monos, "Bad")
})


test_that("validating two non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Bad1", "Bad2")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"

  err <- rlang::catch_cnd(glycan_structure(graph))
  expect_true(grepl("Unknown monosaccharide: Bad1, Bad2", err$parent$message))
  expect_equal(err$parent$monos, c("Bad1", "Bad2"))
})


test_that("validating duplicated non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Bad", "Bad")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"

  err <- rlang::catch_cnd(glycan_structure(graph))
  expect_true(grepl("Unknown monosaccharide: Bad", err$parent$message))
  expect_equal(err$parent$monos, "Bad")
})


test_that("validating bad subtituent", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- c("", "6S", "Bad")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"

  expect_error(glycan_structure(graph))
  err <- rlang::catch_cnd(glycan_structure(graph))
  expect_true(grepl("Unknown substituent: Bad", err$parent$message))
  expect_equal(err$parent$subs, "Bad")
})


patrick::with_parameters_test_that("validating bad linkage", {
    graph <- igraph::make_graph(~ 1-+2, 2-+3)
    igraph::V(graph)$mono <- c("Hex", "Fuc", "Hex")
    igraph::V(graph)$sub <- ""
    igraph::E(graph)$linkage <- bad_linkage
    graph$anomer <- "a1"

    expect_error(glycan_structure(graph))
    err <- rlang::catch_cnd(glycan_structure(graph))
  },
  bad_linkage = c("1-4", "c1-4", "b1", "abc", ""),
  .test_name = bad_linkage
)


test_that("validating NA linkages", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", NA)
  graph$anomer <- "a1"

  expect_error(glycan_structure(graph))
})


test_that("validating mixed generic and concrete monosaccharides", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "GlcNAc", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"

  expect_error(glycan_structure(graph), "Monosaccharides must be either all generic or all concrete")
})


test_that("missing anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"

  expect_error(glycan_structure(graph), "Glycan structure must have a graph attribute 'anomer'")
})


test_that("invalid anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a"

  expect_error(glycan_structure(graph), "Invalid anomer: a")
})




# Tests for glycan_structure vector functions

# Helper function to create a simple glycan structure
create_simple_glycan_graph <- function(mono_names, linkages, anomer = "?1") {
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
  graph
}

# Tests for glycan_structure constructor --------------------------------------------

test_that("glycan_structure creates empty vector by default", {
  sv <- glycan_structure()
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 0)
  expect_length(attr(sv, "structures"), 0)
})

test_that("glycan_structure works with single glycan structure", {
  graph <- o_glycan_core_1()
  sv <- glycan_structure(graph)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 1)
  expect_length(attr(sv, "structures"), 1)
  data <- vctrs::vec_data(sv)
  expect_equal(vctrs::field(data, "iupac")[1], "Gal(b1-3)GalNAc(a1-")
})

test_that("glycan_structure works with multiple different glycan structures", {
  graph1 <- o_glycan_core_1()
  graph2 <- n_glycan_core()
  sv <- glycan_structure(graph1, graph2)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "structures"), 2)
})

test_that("glycan_structure removes duplicates based on IUPAC codes", {
  graph1 <- o_glycan_core_1()
  graph2 <- o_glycan_core_1()  # Same structure
  sv <- glycan_structure(graph1, graph2)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 2)  # Original vector has 2 elements
  expect_length(attr(sv, "structures"), 1)  # But only 1 unique structure
})

test_that("glycan_structure handles structures with different IUPAC but same graph topology", {
  # Create two structures that are topologically the same but have different anomers
  graph1 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4", "a1")
  graph2 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4", "b1")
  
  sv <- glycan_structure(graph1, graph2)
  
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "structures"), 2)  # Different IUPAC codes
})

test_that("glycan_structure validates input", {
  expect_error(glycan_structure("not a graph"), "igraph objects")
  expect_error(glycan_structure(list(1, 2)), "igraph objects")
})

# Tests for as_glycan_structure -------------------------------------------------

test_that("as_glycan_structure creates valid glycan_structure object from igraph", {
  graph <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")
  
  sv <- as_glycan_structure(graph)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 1)
})

test_that("as_glycan_structure works with empty inputs", {
  sv <- glycan_structure()
  result <- as_glycan_structure(sv)
  
  expect_s3_class(result, "glyrepr_structure")
  expect_equal(length(result), 0)
})

# Tests for is_glycan_structure --------------------------------------------------

test_that("is_glycan_structure correctly identifies glycan_structure objects", {
  sv <- glycan_structure()
  expect_true(is_glycan_structure(sv))
  
  expect_false(is_glycan_structure(1))
  expect_false(is_glycan_structure("test"))
  expect_false(is_glycan_structure(list()))
  # Note: individual igraph objects should return FALSE for the new vectorized version
  graph <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")
  expect_false(is_glycan_structure(graph))
})

# Tests for format method -----------------------------------------------------

test_that("format.glyrepr_structure displays correct IUPAC sequences", {
  graph1 <- o_glycan_core_1()
  graph2 <- n_glycan_core()
  sv <- glycan_structure(graph1, graph2)
  
  formatted <- format(sv)
  
  expect_type(formatted, "character")
  expect_length(formatted, 2)
  expect_true("Gal(b1-3)GalNAc(a1-" %in% formatted)
  expect_true("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-" %in% formatted)
})

test_that("format.glyrepr_structure handles empty vector", {
  sv <- glycan_structure()
  formatted <- format(sv)
  
  expect_type(formatted, "character")
  expect_length(formatted, 0)
})

test_that("format.glyrepr_structure handles duplicates correctly", {
  graph1 <- o_glycan_core_1()
  graph2 <- o_glycan_core_1()
  sv <- glycan_structure(graph1, graph2)
  
  formatted <- format(sv)
  
  expect_length(formatted, 2)
  expect_equal(formatted[1], formatted[2])  # Both should show same IUPAC
  expect_equal(formatted[1], "Gal(b1-3)GalNAc(a1-")
})

test_that("truncation works in tibble", {
  sv <- glycan_structure(n_glycan_core(), n_glycan_core(), n_glycan_core())
  tibble <- tibble::tibble(struc = sv, a = 1)
  expect_snapshot(print(tibble, width = 30))
})

# Tests for vctrs methods -----------------------------------------------------

test_that("vec_ptype_abbr.glyrepr_structure returns correct abbreviation", {
  sv <- glycan_structure()
  expect_equal(vctrs::vec_ptype_abbr(sv), "struct")
})

test_that("vec_ptype_full.glyrepr_structure returns correct full type", {
  sv <- glycan_structure()
  expect_equal(vctrs::vec_ptype_full(sv), "glycan_structure")
})

# Tests for print footer ------------------------------------------------------

test_that("obj_print_footer.glyrepr_structure displays unique count", {
  graph1 <- o_glycan_core_1()
  graph2 <- n_glycan_core()
  graph3 <- o_glycan_core_1()  # Duplicate
  sv <- glycan_structure(graph1, graph2, graph3)
  
  output <- capture.output(obj_print_footer.glyrepr_structure(sv))
  
  expect_length(output, 1)
  expect_match(output, "# Unique structures: 2")
})

test_that("obj_print_footer.glyrepr_structure handles empty vector", {
  sv <- glycan_structure()
  
  output <- capture.output(obj_print_footer.glyrepr_structure(sv))
  
  expect_length(output, 1)
  expect_match(output, "# Unique structures: 0")
})

# Tests for get_structure_graphs ----------------------------------------------

test_that("get_structure_graphs extracts individual structures correctly", {
  graph1 <- o_glycan_core_1()
  graph2 <- n_glycan_core()
  sv <- glycan_structure(graph1, graph2)
  
  # Check that we can retrieve the original structures
  extracted_1 <- get_structure_graphs(sv, 1)
  extracted_2 <- get_structure_graphs(sv, 2)
  extracted_all <- get_structure_graphs(sv)
  
  expect_s3_class(extracted_1, "igraph")
  expect_s3_class(extracted_2, "igraph") 
  expect_type(extracted_all, "list")
  expect_length(extracted_all, 2)
})

# Integration tests -----------------------------------------------------------

test_that("glycan_structure preserves glycan structures correctly", {
  graph1 <- o_glycan_core_1()
  graph2 <- n_glycan_core()
  sv <- glycan_structure(graph1, graph2)
  
  # Check that we can retrieve the original structures
  structures <- attr(sv, "structures")
  expect_length(structures, 2)
  
  # Check IUPAC codes are generated correctly
  expected_iupacs <- c(structure_to_iupac(graph1), structure_to_iupac(graph2))
  data <- vctrs::vec_data(sv)
  stored_codes <- vctrs::field(data, "iupac")
  expect_equal(sort(unique(stored_codes)), sort(expected_iupacs))
})

test_that("glycan_structure handles complex branched structures", {
  # Create a more complex branched structure
  complex_graph <- igraph::make_graph(~ 1-+2, 1-+3, 1-+4)
  igraph::V(complex_graph)$mono <- c("Man", "GlcNAc", "Gal", "Fuc")
  igraph::V(complex_graph)$sub <- ""
  igraph::E(complex_graph)$linkage <- c("b1-4", "a1-3", "a1-6")
  complex_graph$anomer <- "a1"
  
  sv <- glycan_structure(complex_graph)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 1)
  expect_length(attr(sv, "structures"), 1)
})

test_that("glycan_structure maintains hash uniqueness property", {
  # Create multiple identical structures
  graphs <- replicate(5, {
    create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")
  }, simplify = FALSE)
  
  sv <- do.call(glycan_structure, graphs)
  
  expect_equal(length(sv), 5)  # 5 elements in vector
  expect_length(attr(sv, "structures"), 1)  # But only 1 unique structure
}) 

# Tests for c() function (vec_ptype2 method) ------------------------------------

test_that("c() combines glycan_structure vectors correctly", {
  # This test ensures the vec_ptype2 method works correctly
  # This was previously failing due to ifelse() with igraph objects
  sv1 <- n_glycan_core()
  sv2 <- o_glycan_core_1()
  
  # This should not error (previously caused rep.igraph error)
  combined <- c(sv1, sv2)
  
  expect_s3_class(combined, "glyrepr_structure")
  expect_equal(length(combined), 2)
  expect_length(attr(combined, "structures"), 2)  # Two unique structures
  
  # Check that both structures are preserved
  formatted <- format(combined)
  expect_true("Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-" %in% formatted)
  expect_true("Gal(b1-3)GalNAc(a1-" %in% formatted)
})

test_that("c() handles duplicate structures across vectors", {
  sv1 <- glycan_structure(o_glycan_core_1(), n_glycan_core())
  sv2 <- glycan_structure(o_glycan_core_1())  # Duplicate of first structure in sv1
  
  combined <- c(sv1, sv2)
  
  expect_s3_class(combined, "glyrepr_structure")
  expect_equal(length(combined), 3)  # Total elements
  expect_length(attr(combined, "structures"), 2)  # Only 2 unique structures
})

test_that("c() works with empty vectors", {
  sv1 <- glycan_structure()
  sv2 <- o_glycan_core_1()
  
  combined1 <- c(sv1, sv2)
  combined2 <- c(sv2, sv1)
  
  expect_s3_class(combined1, "glyrepr_structure")
  expect_s3_class(combined2, "glyrepr_structure")
  expect_equal(length(combined1), 1)
  expect_equal(length(combined2), 1)
  expect_length(attr(combined1, "structures"), 1)
  expect_length(attr(combined2, "structures"), 1)
})

test_that("c() combines multiple structure vectors efficiently", {
  # Test combining multiple vectors with various duplicates
  sv1 <- glycan_structure(o_glycan_core_1())
  sv2 <- glycan_structure(n_glycan_core())
  sv3 <- glycan_structure(o_glycan_core_1(), n_glycan_core())  # Contains both
  
  combined <- c(sv1, sv2, sv3)
  
  expect_s3_class(combined, "glyrepr_structure")
  expect_equal(length(combined), 4)  # 1 + 1 + 2 = 4 total elements
  expect_length(attr(combined, "structures"), 2)  # Only 2 unique structures
  
  # Check all elements are present
  formatted <- format(combined)
  expect_equal(sum(formatted == "Gal(b1-3)GalNAc(a1-"), 2)  # o_glycan_core_1 appears twice
  expect_equal(sum(formatted == "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(?1-"), 2)  # n_glycan_core appears twice
})

test_that("c() preserves structure integrity across combinations", {
  # Use existing example structures
  sv1 <- o_glycan_core_1()
  sv2 <- n_glycan_core()
  
  combined <- c(sv1, sv2)
  
  # Check that we can retrieve the original structures correctly
  extracted_graphs <- get_structure_graphs(combined)
  
  expect_length(extracted_graphs, 2)
  expect_s3_class(extracted_graphs[[1]], "igraph")
  expect_s3_class(extracted_graphs[[2]], "igraph")
  
  # Verify structure properties are preserved using known structures
  # o_glycan_core_1 has GalNAc and Gal
  # n_glycan_core has Man and GlcNAc
  expect_true("GalNAc" %in% igraph::V(extracted_graphs[[1]])$mono)
  expect_true("Gal" %in% igraph::V(extracted_graphs[[1]])$mono)
  expect_true("Man" %in% igraph::V(extracted_graphs[[2]])$mono) 
  expect_true("GlcNAc" %in% igraph::V(extracted_graphs[[2]])$mono)
})

# Tests for vector casting and subsetting functionality -------------------------

test_that("glycan_structure vectors can be subset with necessary structure preservation", {
  sv <- glycan_structure(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  
  # Test various subsetting operations
  subset1 <- sv[1]
  subset2 <- sv[c(1, 3)]
  subset3 <- sv[2:3]
  
  expect_s3_class(subset1, "glyrepr_structure")
  expect_s3_class(subset2, "glyrepr_structure")
  expect_s3_class(subset3, "glyrepr_structure")
  
  expect_equal(length(subset1), 1)
  expect_equal(length(subset2), 2)
  expect_equal(length(subset3), 2)
  
  # Check that unique structure tracking is maintained
  expect_length(attr(sv, "structures"), 2)  # Original has 2 unique
  expect_length(attr(subset1, "structures"), 1)  # Subset preserves 1 unique
  expect_length(attr(subset2, "structures"), 1)  # Subset preserves 1 unique
  expect_length(attr(subset3, "structures"), 2)  # Subset preserves 2 unique
})

test_that("glycan_structure vectors can be repeated", {
  sv <- glycan_structure(o_glycan_core_1())
  
  # Test rep() function which uses vctrs casting methods
  repeated <- rep(sv, 3)
  
  expect_s3_class(repeated, "glyrepr_structure")
  expect_equal(length(repeated), 3)
  expect_length(attr(repeated, "structures"), 1)  # Still only 1 unique structure
  
  formatted <- format(repeated)
  expect_equal(formatted, rep("Gal(b1-3)GalNAc(a1-", 3))
})

test_that("complex vector operations work correctly", {
  # Test more complex vector operations
  sv1 <- glycan_structure(o_glycan_core_1(), n_glycan_core())
  sv2 <- glycan_structure(n_glycan_core(), o_glycan_core_1())
  
  # Combine and then subset
  combined <- c(sv1, sv2)
  reordered <- combined[c(4, 3, 2, 1)]
  
  expect_s3_class(reordered, "glyrepr_structure")
  expect_equal(length(reordered), 4)
  expect_length(attr(reordered, "structures"), 2)
  
  # Check that reordering worked correctly
  formatted_original <- format(combined)
  formatted_reordered <- format(reordered)
  expect_equal(formatted_reordered, rev(formatted_original))
})

test_that("structure vector methods handle edge cases", {
  # Test with single element vector
  single <- glycan_structure(o_glycan_core_1())
  
  # Test combining with itself
  doubled <- c(single, single)
  expect_equal(length(doubled), 2)
  expect_length(attr(doubled, "structures"), 1)
  
  # Test empty + non-empty combinations in different orders
  empty <- glycan_structure()
  combined1 <- c(empty, single, empty)
  combined2 <- c(single, empty, single)
  
  expect_equal(length(combined1), 1)
  expect_equal(length(combined2), 2)
  expect_length(attr(combined1, "structures"), 1)
  expect_length(attr(combined2, "structures"), 1)
})

# Tests for tibble and dplyr operations with structure optimization -----------

test_that("tibble row subsetting optimizes structure storage", {
  skip_if_not_installed("tibble")
  
  # Create test data with duplicates
  sv <- glycan_structure(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv, name = c("A", "B", "C"))
  
  # Original should have 2 unique structures
  expect_length(attr(sv, "structures"), 2)
  
  # Single row subsetting should optimize to 1 unique structure
  subset1 <- df[1, ]
  expect_equal(length(subset1$structure), 1)
  expect_length(attr(subset1$structure, "structures"), 1)
  
  # Multiple row subsetting with same structure should optimize to 1 unique
  subset2 <- df[c(1, 3), ]  # Both point to same structure
  expect_equal(length(subset2$structure), 2)
  expect_length(attr(subset2$structure, "structures"), 1)
  
  # Multiple row subsetting with different structures should keep both
  subset3 <- df[2:3, ]
  expect_equal(length(subset3$structure), 2)
  expect_length(attr(subset3$structure, "structures"), 2)
  
  # Single row with different structure
  subset4 <- df[2, ]
  expect_equal(length(subset4$structure), 1)
  expect_length(attr(subset4$structure, "structures"), 1)
})

test_that("dplyr filter operations optimize structure storage", {
  skip_if_not_installed("dplyr")
  
  # Create test data with duplicates
  sv <- glycan_structure(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv, score = c(10, 20, 30))
  
  # Filter to single row should optimize to 1 unique structure
  filtered1 <- df %>% dplyr::filter(id == 1)
  expect_equal(length(filtered1$structure), 1)
  expect_length(attr(filtered1$structure, "structures"), 1)
  
  # Filter to multiple rows with same structure should optimize
  filtered2 <- df %>% dplyr::filter(id != 2)  # Keeps rows 1 and 3 (same structure)
  expect_equal(length(filtered2$structure), 2)
  expect_length(attr(filtered2$structure, "structures"), 1)
  
  # Filter to multiple rows with different structures should keep both
  filtered3 <- df %>% dplyr::filter(id >= 2)  # Keeps rows 2 and 3 (different structures)
  expect_equal(length(filtered3$structure), 2)
  expect_length(attr(filtered3$structure, "structures"), 2)
  
  # Filter by score
  filtered4 <- df %>% dplyr::filter(score >= 20)
  expect_equal(length(filtered4$structure), 2)
  expect_length(attr(filtered4$structure, "structures"), 2)
})

test_that("dplyr slice operations optimize structure storage", {
  skip_if_not_installed("dplyr")
  
  # Create test data
  sv <- glycan_structure(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv)
  
  # slice() operations
  sliced1 <- df %>% dplyr::slice(1)
  expect_equal(length(sliced1$structure), 1)
  expect_length(attr(sliced1$structure, "structures"), 1)
  
  sliced2 <- df %>% dplyr::slice(c(1, 3))
  expect_equal(length(sliced2$structure), 2)
  expect_length(attr(sliced2$structure, "structures"), 1)
  
  sliced3 <- df %>% dplyr::slice(2:3)
  expect_equal(length(sliced3$structure), 2)
  expect_length(attr(sliced3$structure, "structures"), 2)
  
  # slice_head() and slice_tail()
  head_slice <- df %>% dplyr::slice_head(n = 1)
  expect_equal(length(head_slice$structure), 1)
  expect_length(attr(head_slice$structure, "structures"), 1)
  
  tail_slice <- df %>% dplyr::slice_tail(n = 1)
  expect_equal(length(tail_slice$structure), 1)
  expect_length(attr(tail_slice$structure, "structures"), 1)
})

test_that("dplyr arrange and other operations preserve structure optimization", {
  skip_if_not_installed("dplyr")
  
  # Create test data
  sv <- glycan_structure(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv, score = c(30, 10, 20))
  
  # arrange() should maintain all structures
  arranged <- df %>% dplyr::arrange(score)
  expect_equal(length(arranged$structure), 3)
  expect_length(attr(arranged$structure, "structures"), 2)
  
  # arrange() + slice() should optimize
  arranged_sliced <- df %>% dplyr::arrange(score) %>% dplyr::slice(1)
  expect_equal(length(arranged_sliced$structure), 1)
  expect_length(attr(arranged_sliced$structure, "structures"), 1)
  
  # top_n() operations
  top2 <- df %>% dplyr::top_n(2, score)
  expect_equal(length(top2$structure), 2)
  expect_length(attr(top2$structure, "structures"), 1)  # top_n selects id=1 and id=3, both have same structure
  
  # distinct() operations with duplicated structures
  df_with_dups <- tibble::tibble(
    id = rep(1:3, each = 2),
    structure = rep(sv, each = 2)
  )
  distinct_result <- df_with_dups %>% dplyr::distinct(structure, .keep_all = TRUE)
  expect_equal(length(distinct_result$structure), 2)
  expect_length(attr(distinct_result$structure, "structures"), 2)
})

test_that("complex tibble and dplyr workflows maintain optimization", {
  skip_if_not_installed("dplyr")
  
  # Create more complex test data
  sv <- glycan_structure(
    o_glycan_core_1(), n_glycan_core(), o_glycan_core_1(), 
    n_glycan_core(), o_glycan_core_1()
  )
  df <- tibble::tibble(
    id = 1:5,
    structure = sv,
    type = c("A", "B", "A", "B", "A"),
    score = c(10, 20, 15, 25, 30)
  )
  
  # Original should have 2 unique structures
  expect_length(attr(sv, "structures"), 2)
  
  # Complex workflow: filter + arrange + slice
  result1 <- df %>%
    dplyr::filter(type == "A") %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::slice(1:2)
  
  expect_equal(length(result1$structure), 2)
  expect_length(attr(result1$structure, "structures"), 1)  # All type A have same structure
  
  # Another complex workflow: group operations
  result2 <- df %>%
    dplyr::group_by(type) %>%
    dplyr::slice_max(score, n = 1) %>%
    dplyr::ungroup()
  
  expect_equal(length(result2$structure), 2)
  expect_length(attr(result2$structure, "structures"), 2)  # One from each type
  
  # Workflow with structure column operations
  result3 <- df %>%
    dplyr::filter(score >= 20) %>%
    dplyr::select(structure, score)
  
  expect_equal(length(result3$structure), 3)  # score >= 20 selects 3 rows (id=2,4,5)
  expect_length(attr(result3$structure, "structures"), 2)  # These have 2 different structures
})

test_that("tibble operations preserve structure content integrity", {
  skip_if_not_installed("dplyr")
  
  # Create test data
  sv <- glycan_structure(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv)
  
  # Verify that optimization doesn't affect structure content
  subset_df <- df %>% dplyr::filter(id == 1)
  
  # Extract the structure and verify it's correct
  structure_graph <- get_structure_graphs(subset_df$structure, 1)
  expect_s3_class(structure_graph, "igraph")
  
  # Verify structure content (o_glycan_core_1 has GalNAc and Gal)
  expect_true("GalNAc" %in% igraph::V(structure_graph)$mono)
  expect_true("Gal" %in% igraph::V(structure_graph)$mono)
  
  # Verify IUPAC representation is preserved
  expect_equal(format(subset_df$structure)[1], "Gal(b1-3)GalNAc(a1-")
})

# Tests for vector conversion -------------------------
test_that("converting to character", {
  sv <- glycan_structure(o_glycan_core_1())
  expect_equal(as.character(sv), "Gal(b1-3)GalNAc(a1-")
})
