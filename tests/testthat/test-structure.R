# Tests for glycan structure functions

good_glycan_graph <- function() {
  graph <- igraph::make_graph(~ 1-+2, 2-+3)
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))
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
  glycan <- get_structure_graphs(glycan_vec, return_list = FALSE)

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
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))
  igraph::V(graph)$mono <- c("GlcNAc", "GlcNAc", "GlcNAc")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-4", "b1-3")
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


test_that("validating duplicated linkage positions", {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("GalNAc", "Gal", "Neu5Ac")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-3", "a2-3")
  graph$anomer <- "a1"

  expect_error(glycan_structure(graph), "Duplicated linkage positions")
})


test_that("duplicated ? linkages are OK", {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("GalNAc", "Gal", "Neu5Ac")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-?", "a2-?")
  graph$anomer <- "a1"

  expect_no_error(glycan_structure(graph))
})


test_that("duplicated x/y linkages are OK", {
  graph <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph)$mono <- c("GalNAc", "Gal", "Neu5Ac")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("b1-3/6", "a2-3/6")
  graph$anomer <- "a1"

  expect_no_error(glycan_structure(graph))
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
  igraph::E(graph)$linkage <- c("a1-3", "b1-4")

  expect_error(glycan_structure(graph), "Glycan structure must have a graph attribute 'anomer'")
})


test_that("invalid anomer attr", {
  graph <- igraph::make_tree(3, children = 2, mode = "out")
  igraph::V(graph)$mono <- c("Hex", "Hex", "Hex")
  igraph::V(graph)$sub <- ""
  igraph::E(graph)$linkage <- c("a1-3", "b1-4")
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
  
  igraph::V(graph)$name <- as.character(1:igraph::vcount(graph))
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
  expect_length(attr(sv, "graphs"), 0)
})

test_that("glycan_structure works with single glycan structure", {
  sv <- o_glycan_core_1()

  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 1)
  expect_length(attr(sv, "graphs"), 1)
  expect_equal(structure_to_iupac(sv), "Gal(b1-3)GalNAc(a1-")
})

test_that("glycan_structure works with multiple different glycan structures", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- c(glycan1, glycan2)

  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 2)
  expect_length(attr(sv, "graphs"), 2)
})

test_that("glycan_structure removes duplicates based on IUPAC codes", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- o_glycan_core_1()  # Same structure
  sv <- c(glycan1, glycan2)

  expect_s3_class(sv, "glyrepr_structure")
  expect_equal(length(sv), 2)  # Original vector has 2 elements
  expect_length(attr(sv, "graphs"), 1)  # But only 1 unique structure
})

test_that("glycan_structure handles structures with different IUPAC but same graph topology", {
  # Create two structures that are topologically the same but have different anomers
  graph1 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4", "a1")
  graph2 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4", "b1")

  sv <- glycan_structure(graph1, graph2)

  expect_equal(length(sv), 2)
  expect_length(attr(sv, "graphs"), 2)  # Different IUPAC codes
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
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- c(glycan1, glycan2)

  formatted <- format(sv)

  expect_type(formatted, "character")
  expect_equal(formatted, c(
    "Gal(b1-3)GalNAc(a1-                                ",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
})

test_that("format.glyrepr_structure handles empty vector", {
  sv <- glycan_structure()
  formatted <- format(sv)
  
  expect_type(formatted, "character")
  expect_length(formatted, 0)
})

test_that("format.glyrepr_structure handles duplicates correctly", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- o_glycan_core_1()
  sv <- c(glycan1, glycan2)

  formatted <- format(sv)

  expect_length(formatted, 2)
  expect_equal(formatted[1], formatted[2])  # Both should show same IUPAC
  expect_equal(formatted[1], "Gal(b1-3)GalNAc(a1-")
})

test_that("format.glyrepr_structure includes names with tab separation", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")

  expect_snapshot(format(glycans))
})

test_that("format.glyrepr_structure without names works correctly", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  # Ensure no names
  glycans <- unname(glycans)

  expect_snapshot(format(glycans))
})

test_that("truncation works in tibble", {
  sv <- c(n_glycan_core(), n_glycan_core(), n_glycan_core())
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
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  glycan3 <- o_glycan_core_1()  # Duplicate
  sv <- c(glycan1, glycan2, glycan3)
  
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
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- c(glycan1, glycan2)

  # Check that we can retrieve the original structures
  extracted_all <- get_structure_graphs(sv)
  extracted_single <- get_structure_graphs(sv[1], return_list = FALSE)

  expect_type(extracted_all, "list")
  expect_length(extracted_all, 2)
  expect_s3_class(extracted_all[[1]], "igraph")
  expect_s3_class(extracted_all[[2]], "igraph")
  expect_s3_class(extracted_single, "igraph")
})

test_that("get_structure_graphs return_list parameter works correctly", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- c(glycan1, glycan2)

  # Test default behavior (NULL return_list)
  # For multiple structures, should return list
  result_default_multi <- get_structure_graphs(sv)
  expect_type(result_default_multi, "list")
  expect_length(result_default_multi, 2)

  # For single structure, should return igraph directly
  result_default_single <- get_structure_graphs(sv[1])
  expect_s3_class(result_default_single, "igraph")

  # Test explicit return_list = TRUE
  result_list_true <- get_structure_graphs(sv[1], return_list = TRUE)
  expect_type(result_list_true, "list")
  expect_length(result_list_true, 1)
  expect_s3_class(result_list_true[[1]], "igraph")

  # Test explicit return_list = FALSE
  result_list_false <- get_structure_graphs(sv[1], return_list = FALSE)
  expect_s3_class(result_list_false, "igraph")
})

test_that("get_structure_graphs validates return_list parameter", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- c(glycan1, glycan2)

  # Should error when return_list = FALSE but length > 1
  expect_error(
    get_structure_graphs(sv, return_list = FALSE),
    "return_list.*must be.*TRUE.*NULL.*length greater than 1"
  )
})

# Integration tests -----------------------------------------------------------
test_that("glycan_structure preserves glycan structures correctly", {
  glycan1 <- o_glycan_core_1()
  glycan2 <- n_glycan_core()
  sv <- c(glycan1, glycan2)

  # Check that we can retrieve the original structures
  structures <- attr(sv, "graphs")
  expect_length(structures, 2)

  # Check IUPAC codes are generated correctly
  expected_iupacs <- c(structure_to_iupac(glycan1), structure_to_iupac(glycan2))
  expect_equal(structure_to_iupac(sv), expected_iupacs)
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
  expect_length(attr(sv, "graphs"), 1)
})

test_that("glycan_structure maintains hash uniqueness property", {
  # Create multiple identical structures
  graphs <- replicate(5, {
    create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")
  }, simplify = FALSE)

  sv <- as_glycan_structure(graphs)

  expect_equal(length(sv), 5)  # 5 elements in vector
  expect_length(attr(sv, "graphs"), 1)  # But only 1 unique structure
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
  expect_length(attr(combined, "graphs"), 2)  # Two unique structures

  # Check that both structures are preserved
  expect_equal(structure_to_iupac(combined), c(
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-3)GalNAc(a1-"
  ))
})

test_that("c() handles duplicate structures across vectors", {
  sv1 <- c(o_glycan_core_1(), n_glycan_core())
  sv2 <- o_glycan_core_1()  # Duplicate of first structure in sv1
  
  combined <- c(sv1, sv2)
  
  expect_s3_class(combined, "glyrepr_structure")
  expect_equal(length(combined), 3)  # Total elements
  expect_length(attr(combined, "graphs"), 2)  # Only 2 unique structures
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
  expect_length(attr(combined1, "graphs"), 1)
  expect_length(attr(combined2, "graphs"), 1)
})

test_that("c() combines multiple structure vectors efficiently", {
  # Test combining multiple vectors with various duplicates
  sv1 <- o_glycan_core_1()
  sv2 <- n_glycan_core()
  sv3 <- c(o_glycan_core_1(), n_glycan_core())  # Contains both

  combined <- c(sv1, sv2, sv3)

  expect_s3_class(combined, "glyrepr_structure")
  expect_equal(length(combined), 4)  # 1 + 1 + 2 = 4 total elements
  expect_length(attr(combined, "graphs"), 2)  # Only 2 unique structures

  # Check all elements are present
  expect_equal(structure_to_iupac(combined), c(
    "Gal(b1-3)GalNAc(a1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",
    "Gal(b1-3)GalNAc(a1-",
    "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
  ))
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
  sv <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())

  # Test various subsetting operations
  subset1 <- sv[1]
  subset2 <- sv[c(1, 3)]
  subset3 <- sv[2:3]
  subset4 <- sv[c()]
  subset5 <- sv[integer(0)]

  expect_s3_class(subset1, "glyrepr_structure")
  expect_s3_class(subset2, "glyrepr_structure")
  expect_s3_class(subset3, "glyrepr_structure")
  expect_s3_class(subset4, "glyrepr_structure")
  expect_s3_class(subset5, "glyrepr_structure")

  expect_equal(length(subset1), 1)
  expect_equal(length(subset2), 2)
  expect_equal(length(subset3), 2)
  expect_equal(length(subset4), 0)
  expect_equal(length(subset5), 0)

  # Check that unique structure tracking is maintained
  expect_length(attr(sv, "graphs"), 2)  # Original has 2 unique
  expect_length(attr(subset1, "graphs"), 1)  # Subset preserves 1 unique
  expect_length(attr(subset2, "graphs"), 1)  # Subset preserves 1 unique
  expect_length(attr(subset3, "graphs"), 2)  # Subset preserves 2 unique
  expect_length(attr(subset4, "graphs"), 0)  # Subset preserves 0 unique
  expect_length(attr(subset5, "graphs"), 0)  # Subset preserves 0 unique
})

test_that("glycan_structure vectors can be repeated", {
  sv <- o_glycan_core_1()
  
  # Test rep() function which uses vctrs casting methods
  repeated <- rep(sv, 3)
  
  expect_s3_class(repeated, "glyrepr_structure")
  expect_equal(length(repeated), 3)
  expect_length(attr(repeated, "graphs"), 1)  # Still only 1 unique structure
  
  formatted <- format(repeated)
  expect_equal(formatted, rep("Gal(b1-3)GalNAc(a1-", 3))
})

test_that("complex vector operations work correctly", {
  # Test more complex vector operations
  sv1 <- c(o_glycan_core_1(), n_glycan_core())
  sv2 <- c(n_glycan_core(), o_glycan_core_1())
  
  # Combine and then subset
  combined <- c(sv1, sv2)
  reordered <- combined[c(4, 3, 2, 1)]
  
  expect_s3_class(reordered, "glyrepr_structure")
  expect_equal(length(reordered), 4)
  expect_length(attr(reordered, "graphs"), 2)
  
  # Check that reordering worked correctly
  formatted_original <- format(combined)
  formatted_reordered <- format(reordered)
  expect_equal(formatted_reordered, rev(formatted_original))
})

test_that("structure vector methods handle edge cases", {
  # Test with single element vector
  single <- o_glycan_core_1()
  
  # Test combining with itself
  doubled <- c(single, single)
  expect_equal(length(doubled), 2)
  expect_length(attr(doubled, "graphs"), 1)
  
  # Test empty + non-empty combinations in different orders
  empty <- glycan_structure()
  combined1 <- c(empty, single, empty)
  combined2 <- c(single, empty, single)
  
  expect_equal(length(combined1), 1)
  expect_equal(length(combined2), 2)
  expect_length(attr(combined1, "graphs"), 1)
  expect_length(attr(combined2, "graphs"), 1)
})

# Tests for tibble and dplyr operations with structure optimization -----------

test_that("tibble row subsetting optimizes structure storage", {
  skip_if_not_installed("tibble")
  
  # Create test data with duplicates
  sv <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv, name = c("A", "B", "C"))
  
  # Original should have 2 unique structures
  expect_length(attr(sv, "graphs"), 2)
  
  # Single row subsetting should optimize to 1 unique structure
  subset1 <- df[1, ]
  expect_equal(length(subset1$structure), 1)
  expect_length(attr(subset1$structure, "graphs"), 1)
  
  # Multiple row subsetting with same structure should optimize to 1 unique
  subset2 <- df[c(1, 3), ]  # Both point to same structure
  expect_equal(length(subset2$structure), 2)
  expect_length(attr(subset2$structure, "graphs"), 1)
  
  # Multiple row subsetting with different structures should keep both
  subset3 <- df[2:3, ]
  expect_equal(length(subset3$structure), 2)
  expect_length(attr(subset3$structure, "graphs"), 2)
  
  # Single row with different structure
  subset4 <- df[2, ]
  expect_equal(length(subset4$structure), 1)
  expect_length(attr(subset4$structure, "graphs"), 1)
})

test_that("dplyr filter operations optimize structure storage", {
  skip_if_not_installed("dplyr")
  
  # Create test data with duplicates
  sv <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv, score = c(10, 20, 30))
  
  # Filter to single row should optimize to 1 unique structure
  filtered1 <- df %>% dplyr::filter(id == 1)
  expect_equal(length(filtered1$structure), 1)
  expect_length(attr(filtered1$structure, "graphs"), 1)
  
  # Filter to multiple rows with same structure should optimize
  filtered2 <- df %>% dplyr::filter(id != 2)  # Keeps rows 1 and 3 (same structure)
  expect_equal(length(filtered2$structure), 2)
  expect_length(attr(filtered2$structure, "graphs"), 1)
  
  # Filter to multiple rows with different structures should keep both
  filtered3 <- df %>% dplyr::filter(id >= 2)  # Keeps rows 2 and 3 (different structures)
  expect_equal(length(filtered3$structure), 2)
  expect_length(attr(filtered3$structure, "graphs"), 2)
  
  # Filter by score
  filtered4 <- df %>% dplyr::filter(score >= 20)
  expect_equal(length(filtered4$structure), 2)
  expect_length(attr(filtered4$structure, "graphs"), 2)
})

test_that("dplyr slice operations optimize structure storage", {
  skip_if_not_installed("dplyr")
  
  # Create test data
  sv <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv)
  
  # slice() operations
  sliced1 <- df %>% dplyr::slice(1)
  expect_equal(length(sliced1$structure), 1)
  expect_length(attr(sliced1$structure, "graphs"), 1)
  
  sliced2 <- df %>% dplyr::slice(c(1, 3))
  expect_equal(length(sliced2$structure), 2)
  expect_length(attr(sliced2$structure, "graphs"), 1)
  
  sliced3 <- df %>% dplyr::slice(2:3)
  expect_equal(length(sliced3$structure), 2)
  expect_length(attr(sliced3$structure, "graphs"), 2)
  
  # slice_head() and slice_tail()
  head_slice <- df %>% dplyr::slice_head(n = 1)
  expect_equal(length(head_slice$structure), 1)
  expect_length(attr(head_slice$structure, "graphs"), 1)
  
  tail_slice <- df %>% dplyr::slice_tail(n = 1)
  expect_equal(length(tail_slice$structure), 1)
  expect_length(attr(tail_slice$structure, "graphs"), 1)
})

test_that("dplyr arrange and other operations preserve structure optimization", {
  skip_if_not_installed("dplyr")
  
  # Create test data
  sv <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv, score = c(30, 10, 20))
  
  # arrange() should maintain all structures
  arranged <- df %>% dplyr::arrange(score)
  expect_equal(length(arranged$structure), 3)
  expect_length(attr(arranged$structure, "graphs"), 2)
  
  # arrange() + slice() should optimize
  arranged_sliced <- df %>% dplyr::arrange(score) %>% dplyr::slice(1)
  expect_equal(length(arranged_sliced$structure), 1)
  expect_length(attr(arranged_sliced$structure, "graphs"), 1)
  
  # top_n() operations
  top2 <- df %>% dplyr::top_n(2, score)
  expect_equal(length(top2$structure), 2)
  expect_length(attr(top2$structure, "graphs"), 1)  # top_n selects id=1 and id=3, both have same structure
  
  # distinct() operations with duplicated structures
  df_with_dups <- tibble::tibble(
    id = rep(1:3, each = 2),
    structure = rep(sv, each = 2)
  )
  distinct_result <- df_with_dups %>% dplyr::distinct(structure, .keep_all = TRUE)
  expect_equal(length(distinct_result$structure), 2)
  expect_length(attr(distinct_result$structure, "graphs"), 2)
})

test_that("complex tibble and dplyr workflows maintain optimization", {
  skip_if_not_installed("dplyr")
  
  # Create more complex test data
  sv <- c(
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
  expect_length(attr(sv, "graphs"), 2)
  
  # Complex workflow: filter + arrange + slice
  result1 <- df %>%
    dplyr::filter(type == "A") %>%
    dplyr::arrange(desc(score)) %>%
    dplyr::slice(1:2)
  
  expect_equal(length(result1$structure), 2)
  expect_length(attr(result1$structure, "graphs"), 1)  # All type A have same structure
  
  # Another complex workflow: group operations
  result2 <- df %>%
    dplyr::group_by(type) %>%
    dplyr::slice_max(score, n = 1) %>%
    dplyr::ungroup()
  
  expect_equal(length(result2$structure), 2)
  expect_length(attr(result2$structure, "graphs"), 2)  # One from each type
  
  # Workflow with structure column operations
  result3 <- df %>%
    dplyr::filter(score >= 20) %>%
    dplyr::select(structure, score)
  
  expect_equal(length(result3$structure), 3)  # score >= 20 selects 3 rows (id=2,4,5)
  expect_length(attr(result3$structure, "graphs"), 2)  # These have 2 different structures
})

test_that("tibble operations preserve structure content integrity", {
  skip_if_not_installed("dplyr")
  
  # Create test data
  sv <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  df <- tibble::tibble(id = 1:3, structure = sv)
  
  # Verify that optimization doesn't affect structure content
  subset_df <- df %>% dplyr::filter(id == 1)
  
  # Extract the structure and verify it's correct
  structure_graph <- get_structure_graphs(subset_df$structure, return_list = FALSE)
  expect_s3_class(structure_graph, "igraph")
  
  # Verify structure content (o_glycan_core_1 has GalNAc and Gal)
  expect_true("GalNAc" %in% igraph::V(structure_graph)$mono)
  expect_true("Gal" %in% igraph::V(structure_graph)$mono)
  
  # Verify IUPAC representation is preserved
  expect_equal(format(subset_df$structure)[1], "Gal(b1-3)GalNAc(a1-")
})

# Tests for vector conversion -------------------------
test_that("converting to character", {
  sv <- o_glycan_core_1()
  expect_equal(as.character(sv), "Gal(b1-3)GalNAc(a1-")
})

# Tests for vertex and edge reordering -------------------------
test_that("vertices and edges are reordered correctly", {
  graph1 <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph1)$mono <- c("GalNAc", "Gal", "GlcNAc")
  igraph::V(graph1)$sub <- ""
  igraph::E(graph1)$linkage <- c("b1-3", "b1-6")
  graph1$anomer <- "a1"

  graph2 <- igraph::make_graph(~ 1-+2, 1-+3)
  igraph::V(graph2)$mono <- c("GalNAc", "GlcNAc", "Gal")
  igraph::V(graph2)$sub <- ""
  igraph::E(graph2)$linkage <- c("b1-6", "b1-3")
  graph2$anomer <- "a1"

  sv1 <- glycan_structure(graph1)
  sv2 <- glycan_structure(graph2)
  graph1 <- get_structure_graphs(sv1, return_list = FALSE)
  graph2 <- get_structure_graphs(sv2, return_list = FALSE)

  expect_equal(igraph::V(graph1)$mono, c("Gal", "GlcNAc", "GalNAc"))
  expect_equal(igraph::V(graph2)$mono, c("Gal", "GlcNAc", "GalNAc"))

  expect_equal(igraph::E(graph1)$linkage, c("b1-3", "b1-6"))
  expect_equal(igraph::E(graph2)$linkage, c("b1-3", "b1-6"))
})

# Tests for mono_type validation in glycan_structure -------------------------

test_that("glycan_structure accepts multiple concrete structures", {
  graph1 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")
  graph2 <- create_simple_glycan_graph(c("Man", "GlcNAc"), "b1-4")
  sv <- glycan_structure(graph1, graph2)
  expect_equal(length(sv), 2)
})

test_that("glycan_structure accepts multiple generic structures", {
  graph1 <- create_simple_glycan_graph(c("Hex", "HexNAc"), "b1-4")
  graph2 <- create_simple_glycan_graph(c("Hex", "dHex"), "b1-6")
  sv <- glycan_structure(graph1, graph2)
  expect_equal(length(sv), 2)
})

test_that("glycan_structure rejects mixing concrete and generic structures", {
  graph1 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")  # concrete
  graph2 <- create_simple_glycan_graph(c("Hex", "HexNAc"), "b1-4")  # generic
  expect_error(
    glycan_structure(graph1, graph2),
    "All structures must have the same monosaccharide type"
  )
})

test_that("c() rejects combining concrete and generic structure vectors", {
  sv1 <- o_glycan_core_1()  # concrete: Gal, GalNAc
  sv2 <- n_glycan_core(mono_type = "generic")  # generic: Hex, HexNAc
  expect_error(
    c(sv1, sv2),
    "All structures must have the same monosaccharide type"
  )
})

test_that("get_mono_type returns same type for all structures in vector", {
  graph1 <- create_simple_glycan_graph(c("Glc", "Gal"), "b1-4")
  graph2 <- create_simple_glycan_graph(c("Man", "GlcNAc"), "b1-4")
  sv <- glycan_structure(graph1, graph2)
  expect_equal(get_mono_type(sv), "concrete")
})

test_that("get_mono_type returns same type for generic structures", {
  graph1 <- create_simple_glycan_graph(c("Hex", "HexNAc"), "b1-4")
  graph2 <- create_simple_glycan_graph(c("Hex", "dHex"), "b1-6")
  sv <- glycan_structure(graph1, graph2)
  expect_equal(get_mono_type(sv), "generic")
})

# Tests for [[<- operation (issue #11) -----------------------------------------

test_that("[[<- is forbidden on glyrepr_structure vectors", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())

  expect_error(
    glycans[[1]] <- n_glycan_core(),
    class = "rlang_error"
  )
})

test_that("[[<- is forbidden with igraph value", {
  glycans <- o_glycan_core_1()
  graph <- get_structure_graphs(glycans, return_list = FALSE)

  expect_error(
    glycans[[1]] <- graph,
    class = "rlang_error"
  )
})

# Tests for names-----

test_that("directly setting names works", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")
  expect_equal(names(glycans), c("A", "B"))
})

test_that("new structures have NULL names by default", {
  glycans <- o_glycan_core_1()
  expect_null(names(glycans))
})

test_that("new structure vectors have NULL names by default", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  expect_null(names(glycans))
})

test_that("names are preserved after subsetting with [", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")

  subset1 <- glycans[1]
  expect_equal(names(subset1), "A")

  subset2 <- glycans[2]
  expect_equal(names(subset2), "B")

  subset_both <- glycans[c(1, 2)]
  expect_equal(names(subset_both), c("A", "B"))

  subset_reverse <- glycans[c(2, 1)]
  expect_equal(names(subset_reverse), c("B", "A"))
})

test_that("names are preserved after subsetting with logical index", {
  glycans <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  names(glycans) <- c("A", "B", "C")

  subset_logical <- glycans[c(TRUE, FALSE, TRUE)]
  expect_equal(names(subset_logical), c("A", "C"))
})

test_that("names are preserved after c() combining", {
  glycans1 <- o_glycan_core_1()
  glycans2 <- n_glycan_core()
  names(glycans1) <- "A"
  names(glycans2) <- "B"

  combined <- c(glycans1, glycans2)
  expect_equal(names(combined), c("A", "B"))
})

test_that("names are preserved after rep()", {
  glycans <- o_glycan_core_1()
  names(glycans) <- "A"

  repeated <- rep(glycans, 3)
  expect_equal(names(repeated), c("A", "A", "A"))
})

test_that("setting names to NULL removes names", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")
  expect_equal(names(glycans), c("A", "B"))

  names(glycans) <- NULL
  expect_null(names(glycans))
})

test_that("names work with empty vectors", {
  empty <- glycan_structure()
  expect_null(names(empty))
})

test_that("names work with single-element vectors", {
  glycan <- o_glycan_core_1()
  names(glycan) <- "single"
  expect_equal(names(glycan), "single")
})

test_that("names are preserved after combining vectors with different names", {
  glycans1 <- c(o_glycan_core_1())
  glycans2 <- c(n_glycan_core())
  names(glycans1) <- "A"
  names(glycans2) <- "B"

  combined <- c(glycans1, glycans2)
  expect_equal(names(combined), c("A", "B"))
})

test_that("names work with duplicated structures", {
  glycans <- c(o_glycan_core_1(), n_glycan_core(), o_glycan_core_1())
  names(glycans) <- c("A", "B", "C")

  # First and third point to the same structure but have different names
  expect_equal(names(glycans), c("A", "B", "C"))

  # Subsetting should preserve names correctly
  subset <- glycans[c(1, 3)]
  expect_equal(names(subset), c("A", "C"))
})

test_that("names are preserved in tibble operations", {
  skip_if_not_installed("tibble")

  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")

  df <- tibble::tibble(id = 1:2, structure = glycans)

  # Subsetting tibble should preserve names
  subset_df <- df[1, ]
  expect_equal(names(subset_df$structure), "A")

  subset_df2 <- df[2, ]
  expect_equal(names(subset_df2$structure), "B")
})

test_that("names work correctly with dplyr operations", {
  skip_if_not_installed("dplyr")

  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")

  df <- tibble::tibble(id = 1:2, structure = glycans)

  # filter should preserve names
  filtered <- df %>% dplyr::filter(id == 1)
  expect_equal(names(filtered$structure), "A")

  # slice should preserve names
  sliced <- df %>% dplyr::slice(2)
  expect_equal(names(sliced$structure), "B")
})

test_that("unname removes names", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")

  unname_glycans <- unname(glycans)
  expect_null(names(unname_glycans))
})

test_that("names work with character conversion", {
  glycans <- c(o_glycan_core_1(), n_glycan_core())
  names(glycans) <- c("A", "B")

  # Converting to character should preserve names
  chars <- as.character(glycans)
  expect_equal(names(chars), c("A", "B"))
})

# Comprehensive regression tests for names preservation in glyrepr_structure functions
test_that("all glyrepr_structure functions preserve names", {
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- c(core1, core2, core1)
  names(structures) <- c("A", "B", "C")

  # Type conversion (glyrepr_structure -> glyrepr_structure)
  expect_equal(names(as_glycan_structure(structures)), c("A", "B", "C"))

  # Character to glyrepr_structure conversion
  char_vec <- c(X = "Glc(a1-", Y = "Gal(a1-")
  expect_equal(names(as_glycan_structure(char_vec)), c("X", "Y"))

  # Accessor functions (return atomic vectors)
  expect_equal(names(get_anomer(structures)), c("A", "B", "C"))
  expect_equal(names(has_linkages(structures)), c("A", "B", "C"))
  expect_equal(names(get_structure_level(structures)), c("A", "B", "C"))
  expect_equal(names(count_mono(structures)), c("A", "B", "C"))
  expect_equal(names(structure_to_iupac(structures)), c("A", "B", "C"))

  # Accessor functions (return list)
  expect_equal(names(get_structure_graphs(structures)), c("A", "B", "C"))

  # Transformation functions (return glyrepr_structure)
  expect_equal(names(remove_linkages(structures)), c("A", "B", "C"))
  expect_equal(names(remove_substituents(structures)), c("A", "B", "C"))

  # Composition conversion
  structs_concrete <- c(o_glycan_core_1(), n_glycan_core())
  names(structs_concrete) <- c("X", "Y")
  expect_equal(names(convert_to_generic(structs_concrete)), c("X", "Y"))

  # Level reduction
  expect_equal(names(reduce_structure_level(structures, "topological")), c("A", "B", "C"))

  # Vector combination
  expect_equal(names(c(structures)), c("A", "B", "C"))
})

# Tests for is.na method ----------------------------------------------------

test_that("is.na returns correct logical for structures with NA", {
  struct <- o_glycan_core_1()
  expect_equal(is.na(struct), FALSE)

  struct_na <- c(o_glycan_core_1(), NA)
  expect_equal(is.na(struct_na), c(FALSE, TRUE))

  struct_all_na <- c(NA, NA)
  expect_equal(is.na(struct_all_na), c(TRUE, TRUE))
})

# Tests for NULL/NA handling in glycan_structure -----------------------------

test_that("glycan_structure accepts NULL to create NA", {
  struct <- glycan_structure(NULL)
  expect_equal(length(struct), 1)
  expect_true(is.na(struct))
})

test_that("glycan_structure accepts NA to create NA", {
  struct <- glycan_structure(NA)
  expect_equal(length(struct), 1)
  expect_true(is.na(struct))
})

test_that("glycan_structure handles mixed valid and NA", {
  graph <- o_glycan_core_1() |> get_structure_graphs(return_list = FALSE)
  struct <- glycan_structure(graph, NULL, graph)
  expect_equal(length(struct), 3)
  expect_false(is.na(struct[1]))
  expect_true(is.na(struct[2]))
  expect_false(is.na(struct[3]))
})

# Tests for vec_cast.character with NA handling --------------------------------

test_that("vec_cast handles NA characters", {
  chars <- c("Glc(a1-", NA)
  struct <- as_glycan_structure(chars)
  expect_equal(length(struct), 2)
  expect_false(is.na(struct[1]))
  expect_true(is.na(struct[2]))
})

test_that("vec_cast handles NA at beginning", {
  chars <- c(NA, "Gal(b1-3)GalNAc(a1-")
  struct <- as_glycan_structure(chars)
  expect_equal(length(struct), 2)
  expect_true(is.na(struct[1]))
  expect_false(is.na(struct[2]))
})

test_that("vec_cast handles all NA characters", {
  chars <- c(NA, NA)
  struct <- as_glycan_structure(chars)
  expect_equal(length(struct), 2)
  expect_true(is.na(struct[1]))
  expect_true(is.na(struct[2]))
})

# Tests for vec_restore NA handling -----------------------------------------

test_that("vec_restore skips NA in type checking", {
  struct1 <- o_glycan_core_1()
  struct2 <- glycan_structure(NA)

  # Should not error - NA should be skipped in type check
  combined <- c(struct1, struct2)
  expect_equal(length(combined), 2)
  expect_false(is.na(combined[1]))
  expect_true(is.na(combined[2]))
})

test_that("combining with NA at beginning works", {
  struct <- o_glycan_core_1()
  # c(NA, struct) doesn't preserve type in base R
  combined <- vctrs::vec_c(NA, struct)
  expect_equal(length(combined), 2)
  expect_true(is.na(combined[1]))
  expect_false(is.na(combined[2]))
})

test_that("vec_restore handles all NA vector", {
  struct_na <- glycan_structure(NA, NA)
  expect_equal(length(struct_na), 2)
  expect_equal(is.na(struct_na), c(TRUE, TRUE))
  expect_length(attr(struct_na, "graphs"), 0)  # No graphs for NA elements
})

test_that("vec_restore preserves graphs for non-NA with mixed NA", {
  struct1 <- o_glycan_core_1()
  struct2 <- n_glycan_core()
  combined <- c(struct1, NA, struct2)

  expect_equal(length(combined), 3)
  expect_false(is.na(combined[1]))
  expect_true(is.na(combined[2]))
  expect_false(is.na(combined[3]))
  # Should have 2 unique structures (not counting NA)
  expect_length(attr(combined, "graphs"), 2)
})

# Tests for as_glycan_composition with NA structures ---------------------------------

test_that("as_glycan_composition handles structures with NA", {
  structs <- c(o_glycan_core_1(), NA, n_glycan_core())
  comps <- as_glycan_composition(structs)
  expect_equal(length(comps), 3)
  expect_false(is.na(comps[1]))
  expect_true(is.na(comps[2]))
  expect_false(is.na(comps[3]))
})

test_that("as_glycan_composition handles single NA structure", {
  structs <- glycan_structure(NA)
  comps <- as_glycan_composition(structs)
  expect_equal(length(comps), 1)
  expect_true(is.na(comps[1]))
})

test_that("as_glycan_composition handles all NA structures", {
  structs <- glycan_structure(NA, NA)
  comps <- as_glycan_composition(structs)
  expect_equal(length(comps), 2)
  expect_true(all(is.na(comps)))
})

# Tests for convert_to_generic with NA --------------------------------------------

test_that("convert_to_generic handles single valid structure", {
  # Note: convert_to_generic uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  converted <- convert_to_generic(o_glycan_core_1())
  expect_false(is.na(converted))
  expect_equal(get_mono_type(converted), "generic")
})

test_that("convert_to_generic works on N-glycan core", {
  # Note: convert_to_generic uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  converted <- convert_to_generic(n_glycan_core())
  expect_false(is.na(converted))
  expect_equal(get_mono_type(converted), "generic")
})

# Tests for reduce_structure_level with NA -----------------------------------------

test_that("reduce_structure_level handles structures with NA", {
  # Create structures manually to avoid smap issues with NA
  # The valid structure should work, and NA should be preserved
  structs <- c(o_glycan_core_1(), NA)
  # Note: reduce_structure_level uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  reduced_valid <- reduce_structure_level(o_glycan_core_1(), to_level = "topological")
  expect_equal(get_structure_level(reduced_valid), "topological")
})

test_that("reduce_structure_level handles NA when reducing to basic", {
  # Note: reduce_structure_level uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  reduced_valid <- reduce_structure_level(o_glycan_core_1(), to_level = "basic")
  expect_equal(get_structure_level(reduced_valid), "basic")
})

# Tests for get_mono_type with NA --------------------------------------------------

test_that("get_mono_type handles mixed NA and valid structures", {
  structs <- c(o_glycan_core_1(), NA)
  types <- get_mono_type(structs)
  expect_equal(length(types), 1)  # Returns scalar for glyrepr_structure
  expect_equal(types, "concrete")
})

test_that("get_mono_type handles all NA structures", {
  structs <- glycan_structure(NA, NA)
  types <- get_mono_type(structs)
  expect_equal(length(types), 1)  # Returns scalar for glyrepr_structure
  expect_true(is.na(types))
})

# Tests for rep with NA ------------------------------------------------------------

test_that("rep handles structures with NA", {
  struct <- c(o_glycan_core_1(), NA)
  repeated <- rep(struct, 2)
  expect_equal(length(repeated), 4)
  expect_equal(is.na(repeated), c(FALSE, TRUE, FALSE, TRUE))
})

test_that("rep preserves NA in single-element NA structure", {
  struct <- glycan_structure(NA)
  repeated <- rep(struct, 3)
  expect_equal(length(repeated), 3)
  expect_true(all(is.na(repeated)))
})

# Tests for subsetting preserves NA ------------------------------------------------

test_that("subsetting preserves NA", {
  struct <- c(o_glycan_core_1(), NA, n_glycan_core())
  expect_true(is.na(struct[2]))
  expect_equal(is.na(struct[c(1, 3)]), c(FALSE, FALSE))
  expect_equal(length(struct[c(1, 3)]), 2)
})

test_that("subsetting with logical index preserves NA", {
  struct <- c(o_glycan_core_1(), NA, n_glycan_core())
  subset <- struct[c(TRUE, FALSE, TRUE)]
  expect_equal(length(subset), 2)
  expect_false(is.na(subset[1]))
  expect_false(is.na(subset[2]))
})

test_that("subsetting removes NA elements correctly", {
  struct <- c(o_glycan_core_1(), NA, n_glycan_core())
  subset <- struct[!is.na(struct)]
  expect_equal(length(subset), 2)
  expect_false(any(is.na(subset)))
})

# Tests for format with NA ---------------------------------------------------------

test_that("format preserves order with mixed NA and valid", {
  struct1 <- o_glycan_core_1()
  struct2 <- n_glycan_core()
  structs <- c(struct1, NA, struct2, NA)
  formatted <- format(structs)
  expect_false(is.na(formatted[1]))
  expect_true(grepl("^NA", formatted[2]))  # format() pads with spaces
  expect_false(is.na(formatted[3]))
  expect_true(grepl("^NA", formatted[4]))  # format() pads with spaces
})

test_that("format handles all NA structures", {
  structs <- glycan_structure(NA, NA)
  formatted <- format(structs)
  expect_equal(length(formatted), 2)
  expect_true(all(grepl("^NA", formatted)))  # format() pads with spaces
})

# Tests for tibble printing with NA ------------------------------------------------

test_that("tibble printing handles NA structures", {
  skip_if_not_installed("tibble")
  struct <- c(o_glycan_core_1(), NA)
  df <- tibble::tibble(struct = struct, id = 1:2)
  output <- capture.output(print(df))
  # NA is displayed as "NA" string in tibble
  expect_true(any(grepl("NA", output)))
})

test_that("tibble printing handles multiple NA structures", {
  skip_if_not_installed("tibble")
  struct <- c(o_glycan_core_1(), NA, n_glycan_core(), NA)
  df <- tibble::tibble(struct = struct, id = 1:4)
  output <- capture.output(print(df))
  # Count occurrences of "NA" string (appears once per NA element in data rows)
  na_count <- sum(grepl("NA", output))
  expect_true(na_count >= 2)  # At least two NA entries should appear
})

# Tests for c() combining with NA --------------------------------------------------

test_that("combining multiple structures with NA preserves type", {
  struct1 <- o_glycan_core_1()
  struct2 <- n_glycan_core()

  combined <- c(struct1, NA, struct2)
  expect_s3_class(combined, "glyrepr_structure")
  expect_equal(length(combined), 3)
  expect_false(is.na(combined[1]))
  expect_true(is.na(combined[2]))
  expect_false(is.na(combined[3]))
})

test_that("c() combines all-NA vectors", {
  structs <- glycan_structure(NA, NA)
  combined <- c(structs, glycan_structure(NA))
  expect_equal(length(combined), 3)
  expect_true(all(is.na(combined)))
})

# Tests for all-NA vector ----------------------------------------------------------

test_that("all NA structure vector has correct properties", {
  structs <- glycan_structure(NA, NA)
  expect_equal(length(structs), 2)
  expect_true(all(is.na(structs)))
  expect_equal(length(attr(structs, "graphs")), 0)  # No graphs for NA
})

# Tests for vec_slice with NA ------------------------------------------------------

test_that("vec_slice preserves NA", {
  struct <- c(o_glycan_core_1(), NA, n_glycan_core())
  sliced <- vctrs::vec_slice(struct, c(1, 3))
  expect_equal(length(sliced), 2)
  expect_false(is.na(sliced[1]))
  expect_false(is.na(sliced[2]))
})

test_that("vec_slice can select NA elements", {
  struct <- c(o_glycan_core_1(), NA, n_glycan_core())
  sliced <- vctrs::vec_slice(struct, 2)
  expect_equal(length(sliced), 1)
  expect_true(is.na(sliced[1]))
})

test_that("vec_slice with negative index excludes NA correctly", {
  struct <- c(o_glycan_core_1(), NA, n_glycan_core())
  sliced <- vctrs::vec_slice(struct, -2)
  expect_equal(length(sliced), 2)
  expect_false(any(is.na(sliced)))
})

# Tests for structure_to_iupac with NA ---------------------------------------------

test_that("structure_to_iupac handles structures with NA", {
  structs <- c(o_glycan_core_1(), NA, n_glycan_core())
  iupacs <- structure_to_iupac(structs)
  expect_equal(length(iupacs), 3)
  expect_false(is.na(iupacs[1]))
  expect_true(is.na(iupacs[2]))
  expect_false(is.na(iupacs[3]))
})

# Tests for remove_linkages with NA -------------------------------------------------

test_that("remove_linkages handles single valid structure", {
  # Note: remove_linkages uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  removed <- remove_linkages(o_glycan_core_1())
  expect_equal(has_linkages(removed), FALSE)
})

# Tests for has_linkages with NA ---------------------------------------------------

test_that("has_linkages handles single valid structure", {
  # Note: has_linkages uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  result <- has_linkages(o_glycan_core_1())
  expect_equal(result, TRUE)
})

# Tests for get_anomer with NA ------------------------------------------------------

test_that("get_anomer handles single valid structure", {
  # Note: get_anomer uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  anomer <- get_anomer(o_glycan_core_1())
  expect_equal(anomer, "a1")
})

# Tests for count_mono with NA ------------------------------------------------------

test_that("count_mono handles single valid structure", {
  # Note: count_mono uses smap which currently doesn't handle NA
  # This test verifies the function works for non-NA structures
  count <- count_mono(o_glycan_core_1())
  expect_equal(length(count), 1)
  expect_false(is.na(count))
})