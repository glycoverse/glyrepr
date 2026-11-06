test_that("nested argument permutations stay distinct within one mapping", {
  glycan <- as_glycan_structure("Gal(b1-")
  structures <- rep(glycan, 3)
  values <- list(list(c(1L, 2L, 3L)), list(c(3L, 2L, 1L)), list(c(1L, 2L, 3L)))
  calls <- 0L
  result <- smap2_chr(structures, values, function(graph, data) {
    calls <<- calls + 1L
    paste(data[[1]], collapse = ",")
  })
  expect_identical(result, c("1,2,3", "3,2,1", "1,2,3"))
  expect_identical(calls, 2L)
})

test_that("mapping keys preserve types and missingness", {
  structures <- rep(as_glycan_structure("Gal(b1-"), 3)
  expect_identical(
    smap2_chr(structures, list(1L, 1, "1"), function(graph, value) typeof(value)),
    c("integer", "double", "character")
  )
  expect_identical(
    smap2_lgl(structures[1:2], c(NA_character_, "NA"), function(graph, value) is.na(value)),
    c(TRUE, FALSE)
  )
  names(structures) <- c(NA_character_, "NA", "NA")
  # vctrs normalizes vector names before simap receives them.
  expect_identical(
    unname(simap_chr(structures, function(graph, name) name)),
    names(structures)
  )
})

test_that("argument separators cannot collide across spmap columns", {
  structures <- rep(as_glycan_structure("Gal(b1-"), 2)
  result <- spmap_chr(
    list(structures, c("a|||b", "a"), c("c", "b|||c")),
    function(graph, x, y) paste(x, y, sep = ":")
  )
  expect_identical(result, c("a|||b:c", "a:b|||c"))
})

test_that("structure-returning maps keep distinct nested argument results", {
  structures <- rep(as_glycan_structure("Gal(b1-"), 2)
  values <- list(list(c(1L, 2L)), list(c(2L, 1L)))
  transform <- function(graph, data) {
    igraph::set_vertex_attr(graph, "mono", value = c("Glc", "Man")[data[[1]][1]])
  }
  expect_identical(
    as.character(smap2_structure(structures, values, transform)),
    c("Glc(b1-", "Man(b1-")
  )
  expect_identical(
    as.character(spmap_structure(list(structures, values), transform)),
    c("Glc(b1-", "Man(b1-")
  )
})
