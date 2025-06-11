test_that("smap functions work with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test smap_int with regular function
  result <- smap_int(structures, igraph::vcount)
  expect_equal(length(result), 3)
  expect_type(result, "integer")
  
  # Test smap_chr with regular function
  result_chr <- smap_chr(structures, function(g) g$anomer)
  expect_equal(length(result_chr), 3)
  expect_type(result_chr, "character")
  
  # Test smap_lgl with regular function
  result_lgl <- smap_lgl(structures, function(g) igraph::vcount(g) > 5)
  expect_equal(length(result_lgl), 3)
  expect_type(result_lgl, "logical")
})

test_that("smap functions work with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1() 
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test smap_int with purrr lambda
  result_lambda <- smap_int(structures, ~ igraph::vcount(.x))
  result_regular <- smap_int(structures, igraph::vcount)
  expect_equal(result_lambda, result_regular)
  
  # Test smap_chr with purrr lambda
  result_chr_lambda <- smap_chr(structures, ~ .x$anomer)
  result_chr_regular <- smap_chr(structures, function(g) g$anomer)
  expect_equal(result_chr_lambda, result_chr_regular)
  
  # Test smap_lgl with purrr lambda
  result_lgl_lambda <- smap_lgl(structures, ~ igraph::vcount(.x) > 5)
  result_lgl_regular <- smap_lgl(structures, function(g) igraph::vcount(g) > 5)
  expect_equal(result_lgl_lambda, result_lgl_regular)
})

test_that("smap_unique works with purrr-style lambda functions", {
  # Create test structures with duplicates
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1, core1, core1)
  
  # Test with regular function
  result_regular <- smap_unique(structures, igraph::vcount)
  expect_equal(length(result_regular), 1)  # Only one unique structure
  
  # Test with purrr lambda
  result_lambda <- smap_unique(structures, ~ igraph::vcount(.x))
  expect_equal(result_lambda, result_regular)
})

test_that("smap_structure works with purrr-style lambda functions", {
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
  result_regular <- smap_structure(structures, add_names_regular)
  
  # Test with purrr lambda
  result_lambda <- smap_structure(structures, ~ {
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
  result_lambda <- smap_chr(x, ~ .x$anomer)
  expect_equal(result, result_lambda)
})

test_that("ssome works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()  # smaller structure
  core2 <- n_glycan_core()    # larger structure
  structures <- glycan_structure(core1, core2, core1)
  
  # Test if some structures have more than 5 vertices
  result <- ssome(structures, function(g) igraph::vcount(g) > 5)
  expect_type(result, "logical")
  expect_equal(length(result), 1)
  expect_true(is.logical(result))
  
  # Test if some structures have more than 20 vertices (should be FALSE)
  result_false <- ssome(structures, function(g) igraph::vcount(g) > 20)
  expect_false(result_false)
  
  # Test if some structures have at least 1 vertex (should be TRUE)
  result_true <- ssome(structures, function(g) igraph::vcount(g) >= 1)
  expect_true(result_true)
})

test_that("ssome works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with purrr lambda
  result_lambda <- ssome(structures, ~ igraph::vcount(.x) > 5)
  result_regular <- ssome(structures, function(g) igraph::vcount(g) > 5)
  expect_equal(result_lambda, result_regular)
  
  # Test with another lambda
  result_lambda2 <- ssome(structures, ~ igraph::vcount(.x) >= 1)
  expect_true(result_lambda2)
})

test_that("severy works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test if all structures have at least 1 vertex (should be TRUE)
  result_true <- severy(structures, function(g) igraph::vcount(g) >= 1)
  expect_true(result_true)
  
  # Test if all structures have more than 10 vertices (should be FALSE)
  result_false <- severy(structures, function(g) igraph::vcount(g) > 10)
  expect_false(result_false)
  
  # Test if all structures have at least 3 vertices
  result <- severy(structures, function(g) igraph::vcount(g) >= 3)
  expect_type(result, "logical")
  expect_equal(length(result), 1)
})

test_that("severy works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with purrr lambda
  result_lambda <- severy(structures, ~ igraph::vcount(.x) >= 1)
  result_regular <- severy(structures, function(g) igraph::vcount(g) >= 1)
  expect_equal(result_lambda, result_regular)
  
  # Test with another lambda
  result_lambda2 <- severy(structures, ~ igraph::vcount(.x) >= 3)
  expect_type(result_lambda2, "logical")
})

test_that("snone works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test if no structures have more than 50 vertices (should be TRUE)
  result_true <- snone(structures, function(g) igraph::vcount(g) > 50)
  expect_true(result_true)
  
  # Test if no structures have at least 1 vertex (should be FALSE)
  result_false <- snone(structures, function(g) igraph::vcount(g) >= 1)
  expect_false(result_false)
  
  # Test return type
  result <- snone(structures, function(g) igraph::vcount(g) > 20)
  expect_type(result, "logical")
  expect_equal(length(result), 1)
})

test_that("snone works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with purrr lambda
  result_lambda <- snone(structures, ~ igraph::vcount(.x) > 50)
  result_regular <- snone(structures, function(g) igraph::vcount(g) > 50)
  expect_equal(result_lambda, result_regular)
  
  # Test with another lambda
  result_lambda2 <- snone(structures, ~ igraph::vcount(.x) >= 1)
  expect_false(result_lambda2)
})

test_that("structure predicate functions work with duplicate structures efficiently", {
  # Create structures with many duplicates to test efficiency
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1, core1, core1, core1, core1)  # same structure 5 times
  
  # These should only evaluate the predicate once for the unique structure
  result_some <- ssome(structures, function(g) igraph::vcount(g) > 3)
  result_every <- severy(structures, function(g) igraph::vcount(g) > 3)
  result_none <- snone(structures, function(g) igraph::vcount(g) > 50)
  
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
  expect_type(ssome(single, ~ igraph::vcount(.x) > 0), "logical")
  expect_type(severy(single, ~ igraph::vcount(.x) > 0), "logical")
  expect_type(snone(single, ~ igraph::vcount(.x) > 50), "logical")
})

# Tests for smap2 functions -----------------------------------------

test_that("smap2 functions work with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  weights <- c(1.0, 2.0, 1.5)
  
  # Test smap2_dbl with regular function
  result <- smap2_dbl(structures, weights, function(g, w) igraph::vcount(g) * w)
  expect_equal(length(result), 3)
  expect_type(result, "double")
  
  # Test smap2_int with regular function
  n_vertices <- c(1, 2, 1)
  result_int <- smap2_int(structures, n_vertices, function(g, n) igraph::vcount(g) + n)
  expect_equal(length(result_int), 3)
  expect_type(result_int, "integer")
  
  # Test smap2_lgl with regular function
  thresholds <- c(5, 6, 5)
  result_lgl <- smap2_lgl(structures, thresholds, function(g, t) igraph::vcount(g) > t)
  expect_equal(length(result_lgl), 3)
  expect_type(result_lgl, "logical")
})

test_that("smap2 functions work with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  weights <- c(1.0, 2.0, 1.5)
  
  # Test smap2_dbl with purrr lambda
  result_lambda <- smap2_dbl(structures, weights, ~ igraph::vcount(.x) * .y)
  result_regular <- smap2_dbl(structures, weights, function(g, w) igraph::vcount(g) * w)
  expect_equal(result_lambda, result_regular)
  
  # Test smap2_lgl with purrr lambda
  thresholds <- c(5, 6, 5)
  result_lgl_lambda <- smap2_lgl(structures, thresholds, ~ igraph::vcount(.x) > .y)
  result_lgl_regular <- smap2_lgl(structures, thresholds, function(g, t) igraph::vcount(g) > t)
  expect_equal(result_lgl_lambda, result_lgl_regular)
})

test_that("smap2 functions work with recycling", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  
  # Test with single weight (should be recycled)
  single_weight <- 2.5
  result <- smap2_dbl(structures, single_weight, ~ igraph::vcount(.x) * .y)
  expect_equal(length(result), 3)
  expect_type(result, "double")
  
  # All results should use the same weight
  expected <- smap_dbl(structures, ~ igraph::vcount(.x) * 2.5)
  expect_equal(result, expected)
})

test_that("smap2_structure works with regular functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2)  # different structures
  
  # Function that adds a graph attribute (doesn't change structure topology)
  add_attribute <- function(g, attr_val) {
    igraph::set_graph_attr(g, "test_attr", attr_val)
  }
  
  attr_values <- c("value1", "value2")
  result <- smap2_structure(structures, attr_values, add_attribute)
  
  expect_s3_class(result, "glyrepr_structure")
  expect_equal(length(result), 2)
  
  # Check that the attributes were added correctly
  result_graphs <- get_structure_graphs(result)
  expect_equal(igraph::graph_attr(result_graphs[[1]], "test_attr"), "value1")
  expect_equal(igraph::graph_attr(result_graphs[[2]], "test_attr"), "value2")
})

test_that("smap2_structure works with purrr-style lambda functions", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2)  # different structures
  
  # Test with purrr lambda that adds vertex attribute based on second argument
  labels <- c("label1", "label2")
  result <- smap2_structure(structures, labels, ~ {
    igraph::set_graph_attr(.x, "custom_label", .y)
  })
  
  expect_s3_class(result, "glyrepr_structure")
  expect_equal(length(result), 2)
  
  # Check that the custom attributes were added
  result_graphs <- get_structure_graphs(result)
  expect_equal(igraph::graph_attr(result_graphs[[1]], "custom_label"), "label1")
  expect_equal(igraph::graph_attr(result_graphs[[2]], "custom_label"), "label2")
})

test_that("smap2 functions handle duplicate structures efficiently", {
  # Create structures with duplicates to test efficiency
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1, core1, core1)  # same structure 3 times
  weights <- c(1.0, 2.0, 1.0)  # first and third are same
  
  # This should only compute twice: once for (core1, 1.0) and once for (core1, 2.0)
  result <- smap2_dbl(structures, weights, function(g, w) igraph::vcount(g) * w)
  
  expect_equal(length(result), 3)
  expect_type(result, "double")
  
  # First and third results should be equal since they have same structure and weight
  expect_equal(result[1], result[3])
  # Second result should be different
  expect_true(result[2] != result[1])
})

test_that("smap2 functions validate inputs", {
  core1 <- o_glycan_core_1()
  structures <- glycan_structure(core1)
  
  # Test that first argument must be glycan_structure
  expect_error(smap2_dbl(c(1, 2, 3), c(1, 2, 3), ~ .x * .y), "glycan_structure")
  expect_error(smap2_structure("not_a_structure", c(1), ~ .x), "glycan_structure")
  
  # Test that smap2_structure validates return type
  expect_error(
    smap2_structure(structures, 1, ~ "not_an_igraph"), 
    "igraph object"
  )
})

test_that("smap2 functions handle edge cases", {
  # Test with single structure
  single <- glycan_structure(o_glycan_core_1())
  result <- smap2_dbl(single, 3.0, ~ igraph::vcount(.x) * .y)
  expect_equal(length(result), 1)
  expect_type(result, "double")
  
  # Test with empty structures
  empty <- glycan_structure()
  empty_result <- smap2_dbl(empty, numeric(0), ~ igraph::vcount(.x) * .y)
  expect_equal(length(empty_result), 0)
  expect_type(empty_result, "double")
})

test_that("smap2 basic functionality works", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  values <- c("a", "b", "c")
  
  # Test smap2 (returns list)
  result <- smap2(structures, values, function(g, v) {
    list(vcount = igraph::vcount(g), value = v)
  })
  
  expect_type(result, "list")
  expect_equal(length(result), 3)
  expect_equal(result[[1]]$value, "a")
  expect_equal(result[[2]]$value, "b")
  expect_equal(result[[3]]$value, "c")
})

test_that("smap2_chr works correctly", {
  # Create test structures
  core1 <- o_glycan_core_1()
  core2 <- n_glycan_core()
  structures <- glycan_structure(core1, core2, core1)
  prefixes <- c("prefix1_", "prefix2_", "prefix3_")
  
  # Test smap2_chr
  result <- smap2_chr(structures, prefixes, function(g, p) {
    paste0(p, igraph::vcount(g))
  })
  
  expect_type(result, "character")
  expect_equal(length(result), 3)
  expect_true(all(grepl("^prefix[123]_\\d+$", result)))
}) 