test_that("structure_map functions work with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test structure_map_int with regular function
  result <- structure_map_int(structures, igraph::vcount)
  expect_equal(length(result), 3)
  expect_type(result, "integer")
  
  # Test structure_map_chr with regular function
  result_chr <- structure_map_chr(structures, function(g) g$anomer)
  expect_equal(length(result_chr), 3)
  expect_type(result_chr, "character")
  
  # Test structure_map_lgl with regular function
  result_lgl <- structure_map_lgl(structures, function(g) igraph::vcount(g) > 5)
  expect_equal(length(result_lgl), 3)
  expect_type(result_lgl, "logical")
})

test_that("structure_map functions work with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1() 
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test structure_map_int with purrr lambda
  result_lambda <- structure_map_int(structures, ~ igraph::vcount(.x))
  result_regular <- structure_map_int(structures, igraph::vcount)
  expect_equal(result_lambda, result_regular)
  
  # Test structure_map_chr with purrr lambda
  result_chr_lambda <- structure_map_chr(structures, ~ .x$anomer)
  result_chr_regular <- structure_map_chr(structures, function(g) g$anomer)
  expect_equal(result_chr_lambda, result_chr_regular)
  
  # Test structure_map_lgl with purrr lambda
  result_lgl_lambda <- structure_map_lgl(structures, ~ igraph::vcount(.x) > 5)
  result_lgl_regular <- structure_map_lgl(structures, function(g) igraph::vcount(g) > 5)
  expect_equal(result_lgl_lambda, result_lgl_regular)
})

test_that("structure_map_unique works with purrr-style lambda functions", {
  # Create test structures with duplicates
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1, core1, core1)
  
  # Test with regular function
  result_regular <- structure_map_unique(structures, igraph::vcount)
  expect_equal(length(result_regular), 1)  # Only one unique structure
  
  # Test with purrr lambda
  result_lambda <- structure_map_unique(structures, ~ igraph::vcount(.x))
  expect_equal(result_lambda, result_regular)
})

test_that("structure_map_structure works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1, core1)
  
  # Function that adds vertex names if not present
  add_names_regular <- function(g) {
    if (!("name" %in% igraph::vertex_attr_names(g))) {
      igraph::set_vertex_attr(g, "name", value = paste0("v", seq_len(igraph::vcount(g))))
    } else {
      g
    }
  }
  
  # Test with regular function
  result_regular <- structure_map_structure(structures, add_names_regular)
  
  # Test with purrr lambda
  result_lambda <- structure_map_structure(structures, ~ {
    if (!("name" %in% igraph::vertex_attr_names(.x))) {
      igraph::set_vertex_attr(.x, "name", value = paste0("v", seq_len(igraph::vcount(.x))))
    } else {
      .x
    }
  })
  
  # Both results should be equivalent
  expect_s3_class(result_lambda, "glyrepr_structure")
  expect_equal(length(result_lambda), length(result_regular))
})

test_that("get_anomer function works correctly", {
  # Test the specific case from the user's example
  x <- n_glycan_core()
  result <- get_anomer(x)
  expect_type(result, "character")
  expect_equal(length(result), 1)
  expect_equal(result, "?1")
  
  # Test that it's equivalent to purrr lambda style
  result_lambda <- structure_map_chr(x, ~ .x$anomer)
  expect_equal(result, result_lambda)
}) 