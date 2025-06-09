# Tests for glycan structure functions

good_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Glc", "Gal", "Glc")
  igraph::V(graph)$sub <- c("", "", "")
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE
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


# Tests for validate_glycan_structure --------------------------------------------

test_that("validate_glycan_structure accepts valid graphs", {
  graph <- good_glycan_graph()
  expect_no_error(validate_glycan_structure(graph))
})


test_that("validate_glycan_structure rejects invalid graphs", {
  # Test undirected graph
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("Glc", "Gal")
  igraph::V(graph)$sub <- c("", "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
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
  graph$alditol <- FALSE
  glycan <- glycan_structure(graph)
  expect_s3_class(glycan, c("glyrepr_structure"))
})


test_that("validating undirected graphs", {
  graph <- igraph::make_graph(~ 1--2)
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must be directed")
})


test_that("validating an in tree", {
  graph <- igraph::make_tree(3, children = 2, mode = "in")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must be an out tree")
})


test_that("validating graph without monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must have a vertex attribute 'mono'")
})


test_that("validating graph without substituent attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must have a vertex attribute 'sub'")
})


test_that("validating graph with NA in monosaccharide attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", NA, "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must have no NA in vertex attribute 'mono'")
})


test_that("validating graph with NA in substitude attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- c("", NA, "")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must have no NA in vertex attribute 'sub'")
})


patrick::with_parameters_test_that("valid substituents", {
  skip_on_old_win()
  graph <- igraph::make_empty_graph(n = 1)
  igraph::V(graph)$mono <- "GlcNAc"
  igraph::V(graph)$sub <- sub
  igraph::E(graph)$linkage <- character(0)
  graph$anomer <- "b1"
  graph$alditol <- FALSE
  expect_no_error(glycan_structure(graph))
}, sub = c("6S", "9Ac", "2P", "?S"))


test_that("validating graph without linkage attribute", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  expect_error(glycan_structure(graph), "Glycan structure must have an edge attribute 'linkage'")
})


test_that("validating one non-existing monosaccharide", {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$mono <- c("Hex", "Fuc", "Bad")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-4")
  graph$anomer <- "a1"
  graph$alditol <- FALSE

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
  graph$alditol <- FALSE

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
  graph$alditol <- FALSE

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
  graph$alditol <- FALSE

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
    graph$alditol <- FALSE

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
  graph$alditol <- FALSE

  expect_error(glycan_structure(graph))
})


test_that("validating mixed generic and concrete monosaccharides", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "GlcNAc", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE

  expect_error(glycan_structure(graph), "Monosaccharides must be either all generic or all concrete")
})


test_that("missing anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$alditol <- FALSE

  expect_error(glycan_structure(graph), "Glycan structure must have a graph attribute 'anomer'")
})


test_that("invalid anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a"
  graph$alditol <- FALSE

  expect_error(glycan_structure(graph), "Invalid anomer: a")
})


test_that("missing alditol attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"

  expect_error(glycan_structure(graph), "Glycan structure must have a graph attribute 'alditol'")
})


test_that("invalid alditol attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- "True"

  expect_error(glycan_structure(graph), "Glycan structure attribute 'alditol' must be logical")
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
  graph$alditol <- FALSE
  graph
}

# Tests for glycan_structure constructor --------------------------------------------

test_that("glycan_structure creates empty vector by default", {
  sv <- glycan_structure()
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 0)
  expect_length(attr(sv, "iupacs"), 0)
  expect_length(attr(sv, "structures"), 0)
})

test_that("glycan_structure works with single glycan structure", {
  graph <- o_glycan_core_1()
  sv <- glycan_structure(graph)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 1)
  expect_length(attr(sv, "iupacs"), 1)
  expect_length(attr(sv, "structures"), 1)
  expect_equal(attr(sv, "iupacs")[[1]], "Gal(b1-3)GalNAc(a1-")
})

test_that("glycan_structure works with multiple different glycan structures", {
  graph1 <- o_glycan_core_1()
  graph2 <- n_glycan_core()
  sv <- glycan_structure(graph1, graph2)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "iupacs"), 2)
  expect_length(attr(sv, "structures"), 2)
})

test_that("glycan_structure removes duplicates based on IUPAC codes", {
  graph1 <- o_glycan_core_1()
  graph2 <- o_glycan_core_1()  # Same structure
  sv <- glycan_structure(graph1, graph2)
  
  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 2)  # Original vector has 2 elements
  expect_length(attr(sv, "iupacs"), 1)  # But only 1 unique structure
  expect_length(attr(sv, "structures"), 1)
})

test_that("glycan_structure handles structures with different IUPAC but same graph topology", {
  # Create two structures that are topologically the same but have different anomers
  graph1 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4", "a1")
  graph2 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4", "b1")
  
  sv <- glycan_structure(graph1, graph2)
  
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "iupacs"), 2)  # Different IUPAC codes
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

# Tests for vctrs methods -----------------------------------------------------

test_that("vec_ptype_abbr.glyrepr_structure returns correct abbreviation", {
  sv <- glycan_structure()
  expect_equal(vctrs::vec_ptype_abbr(sv), "structure")
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
  stored_iupacs <- attr(sv, "iupacs")
  expect_equal(sort(unname(stored_iupacs)), sort(expected_iupacs))
})

test_that("glycan_structure handles complex branched structures", {
  # Create a more complex branched structure
  complex_graph <- igraph::make_graph(~ 1-+2, 1-+3, 1-+4)
  igraph::V(complex_graph)$mono <- c("Man", "GlcNAc", "Gal", "Fuc")
  igraph::V(complex_graph)$sub <- ""
  igraph::E(complex_graph)$linkage <- c("b1-4", "a1-3", "a1-6")
  complex_graph$anomer <- "a1"
  complex_graph$alditol <- FALSE
  
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
  expect_length(attr(sv, "iupacs"), 1)
}) 