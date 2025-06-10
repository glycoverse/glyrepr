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

test_that("structure_some works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()  # smaller structure
  core2 <- n_glycan_core()    # larger structure
  structures <- glycan_structure(core1, core2, core1)
  
  # Test if some structures have more than 5 vertices
  result <- structure_some(structures, function(g) igraph::vcount(g) > 5)
  expect_type(result, "logical")
  expect_equal(length(result), 1)
  expect_true(is.logical(result))
  
  # Test if some structures have more than 20 vertices (should be FALSE)
  result_false <- structure_some(structures, function(g) igraph::vcount(g) > 20)
  expect_false(result_false)
  
  # Test if some structures have at least 1 vertex (should be TRUE)
  result_true <- structure_some(structures, function(g) igraph::vcount(g) >= 1)
  expect_true(result_true)
})

test_that("structure_some works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with purrr lambda
  result_lambda <- structure_some(structures, ~ igraph::vcount(.x) > 5)
  result_regular <- structure_some(structures, function(g) igraph::vcount(g) > 5)
  expect_equal(result_lambda, result_regular)
  
  # Test with another lambda
  result_lambda2 <- structure_some(structures, ~ igraph::vcount(.x) >= 1)
  expect_true(result_lambda2)
})

test_that("structure_every works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test if all structures have at least 1 vertex (should be TRUE)
  result_true <- structure_every(structures, function(g) igraph::vcount(g) >= 1)
  expect_true(result_true)
  
  # Test if all structures have more than 10 vertices (should be FALSE)
  result_false <- structure_every(structures, function(g) igraph::vcount(g) > 10)
  expect_false(result_false)
  
  # Test if all structures have at least 3 vertices
  result <- structure_every(structures, function(g) igraph::vcount(g) >= 3)
  expect_type(result, "logical")
  expect_equal(length(result), 1)
})

test_that("structure_every works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with purrr lambda
  result_lambda <- structure_every(structures, ~ igraph::vcount(.x) >= 1)
  result_regular <- structure_every(structures, function(g) igraph::vcount(g) >= 1)
  expect_equal(result_lambda, result_regular)
  
  # Test with another lambda
  result_lambda2 <- structure_every(structures, ~ igraph::vcount(.x) >= 3)
  expect_type(result_lambda2, "logical")
})

test_that("structure_none works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test if no structures have more than 50 vertices (should be TRUE)
  result_true <- structure_none(structures, function(g) igraph::vcount(g) > 50)
  expect_true(result_true)
  
  # Test if no structures have at least 1 vertex (should be FALSE)
  result_false <- structure_none(structures, function(g) igraph::vcount(g) >= 1)
  expect_false(result_false)
  
  # Test return type
  result <- structure_none(structures, function(g) igraph::vcount(g) > 20)
  expect_type(result, "logical")
  expect_equal(length(result), 1)
})

test_that("structure_none works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with purrr lambda
  result_lambda <- structure_none(structures, ~ igraph::vcount(.x) > 50)
  result_regular <- structure_none(structures, function(g) igraph::vcount(g) > 50)
  expect_equal(result_lambda, result_regular)
  
  # Test with another lambda
  result_lambda2 <- structure_none(structures, ~ igraph::vcount(.x) >= 1)
  expect_false(result_lambda2)
})

test_that("structure predicate functions work with duplicate structures efficiently", {
  # Create structures with many duplicates to test efficiency
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1, core1, core1, core1, core1)  # same structure 5 times
  
  # These should only evaluate the predicate once for the unique structure
  result_some <- structure_some(structures, function(g) igraph::vcount(g) > 3)
  result_every <- structure_every(structures, function(g) igraph::vcount(g) > 3)
  result_none <- structure_none(structures, function(g) igraph::vcount(g) > 50)
  
  expect_type(result_some, "logical")
  expect_type(result_every, "logical")
  expect_type(result_none, "logical")
  
  expect_equal(length(result_some), 1)
  expect_equal(length(result_every), 1)
  expect_equal(length(result_none), 1)
})

test_that("structure predicate functions handle edge cases", {
  # Test with single structure
  single <- glycan_structure(o_glycan_core_1())
  expect_type(structure_some(single, ~ igraph::vcount(.x) > 0), "logical")
  expect_type(structure_every(single, ~ igraph::vcount(.x) > 0), "logical")
  expect_type(structure_none(single, ~ igraph::vcount(.x) > 50), "logical")
}) 