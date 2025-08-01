# Tests for has_linkages() function
test_that("`has_linkages()` works for glycans with intact linkages", {
  glycan <- o_glycan_core_1(linkage = TRUE)
  expect_true(has_linkages(glycan))
})


test_that("`has_linkages()` works for glycan with partial linkages", {
  glycan_vec <- o_glycan_core_2(linkage = TRUE)
  glycan <- get_structure_graphs(glycan_vec, 1)  # Extract the igraph
  glycan <- igraph::set_edge_attr(glycan, "linkage", value = c("??-?", "b1-6"))
  
  # Create a new vectorized structure with the modified graph
  modified_glycan_vec <- glycan_structure(glycan)
  expect_true(has_linkages(modified_glycan_vec))
})


test_that("`has_linkages()` works for glycans without linkages", {
  glycan <- o_glycan_core_1(linkage = FALSE)
  expect_false(has_linkages(glycan))
})


# Tests for possible_linkages() function
test_that("a?-2", {
  result <- possible_linkages("a?-2")
  expected <- c("a1-2", "a2-2")
  expect_identical(result, expected)
})


test_that("a1-?", {
  result <- possible_linkages("a1-?")
  expected <- paste0("a1-", 1:9)
  expect_identical(result, expected)
})


test_that("?1-3", {
  result <- possible_linkages("?1-3")
  expected <- c("a1-3", "b1-3")
  expect_identical(result, expected)
})


test_that("b?-?", {
  result <- possible_linkages("b?-?")
  df <- expand.grid(1:2, 1:9)
  expected <- apply(df, 1, function(x) paste0("b", x[1], "-", x[2]))
  expect_identical(result, expected)
})


test_that("a2-3|6", {
  result <- possible_linkages("a2-3/6")
  expected <- c("a2-3", "a2-6")
  expect_identical(result, expected)
})


test_that("a1-3", {
  expect_identical(possible_linkages("a1-3"), "a1-3")
})


test_that("wrong format", {
  expect_error(possible_linkages("a1-3-4"), "Invalid linkage format")
  expect_error(possible_linkages(""), "Invalid linkage format")
})


test_that("wrong input type", {
  expect_error(possible_linkages(1))
  expect_error(possible_linkages(NA))
  expect_error(possible_linkages(NULL))
})


test_that("multiple linkages", {
  expect_error(possible_linkages(c("a1-3", "b1-4")))
})


test_that("custom ranges", {
  result <- possible_linkages("??-?", anomer_range = "a", pos1_range = 2, pos2_range = c(3, 6))
  expected <- c("a2-3", "a2-6")
  expect_identical(result, expected)
})


test_that("custom ranges wrong types", {
  expect_error(possible_linkages("??-?", anomer_range = 1))
  expect_error(possible_linkages("??-?", pos1_range = "a"))
  expect_error(possible_linkages("??-?", pos2_range = "a"))
})


test_that("0 range", {
  expect_equal(possible_linkages("??-?", pos1_range = integer(0)), character(0))
})


test_that("include unknown", {
  result <- possible_linkages("?1-6", include_unknown = TRUE)
  expected <- c("a1-6", "b1-6", "?1-6")
  expect_identical(result, expected)
})


# Tests for remove_linkages() function
test_that("works for glycans", {
  glycan <- o_glycan_core_2(linkage = TRUE)
  glycan <- remove_linkages(glycan)
  expect_false(has_linkages(glycan))
})
