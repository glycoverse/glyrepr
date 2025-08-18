# Tests for substituent functions

test_that("available_substituents returns expected values", {
  subs <- available_substituents()
  expect_type(subs, "character")
  expect_true(length(subs) > 0)
  expect_true("Me" %in% subs)
  expect_true("Ac" %in% subs)
  expect_true("S" %in% subs)
})

test_that("normalize_substituents handles empty strings", {
  expect_equal(normalize_substituents(""), "")
})

test_that("normalize_substituents handles single substituents", {
  expect_equal(normalize_substituents("6S"), "6S")
  expect_equal(normalize_substituents("3Me"), "3Me")
  expect_equal(normalize_substituents("?Ac"), "?Ac")
})

test_that("normalize_substituents sorts multiple substituents by position", {
  expect_equal(normalize_substituents("4Ac,3Me"), "3Me,4Ac")
  expect_equal(normalize_substituents("6P,2S,4Ac"), "2S,4Ac,6P")
  expect_equal(normalize_substituents("9Ac,1Me,5S"), "1Me,5S,9Ac")
})

test_that("normalize_substituents handles ? positions correctly", {
  expect_equal(normalize_substituents("?S,3Me"), "3Me,?S")
  expect_equal(normalize_substituents("?Ac,?Me"), "?Ac,?Me")  # Original order preserved for same position
  expect_equal(normalize_substituents("4Ac,?S,2Me"), "2Me,4Ac,?S")
})

test_that("normalize_substituents removes empty components", {
  expect_equal(normalize_substituents("3Me,,4Ac"), "3Me,4Ac")
  expect_equal(normalize_substituents(",6S,"), "6S")
})

test_that("normalize_substituents input validation", {
  expect_error(normalize_substituents(c("3Me", "4Ac")))
  expect_error(normalize_substituents(123))
})

test_that("remove_substituents works on glycan structures", {
  # Create a test glycan with substituents
  graph <- igraph::make_graph(~ 1-+2)
  igraph::V(graph)$mono <- c("Glc", "Gal")
  igraph::V(graph)$sub <- c("3Me,4Ac", "6S")
  igraph::E(graph)$linkage <- "b1-4"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  
  glycan <- glycan_structure(graph)
  
  # Remove substituents
  clean_glycan <- remove_substituents(glycan)
  
  # Check that substituents are removed
  clean_graph <- get_structure_graphs(clean_glycan, return_list = FALSE)
  expect_equal(igraph::V(clean_graph)$sub, c("", ""))
  expect_equal(igraph::V(clean_graph)$mono, c("Gal", "Glc"))  # Monos should be unchanged
})

test_that("remove_substituents input validation", {
  expect_error(remove_substituents("not_a_glycan"))
  expect_error(remove_substituents(123))
})
