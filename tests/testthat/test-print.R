test_that("print works for NE glycan graphs", {
  skip_on_old_win()
  x <- new_ne_glycan_graph(n_glycan_core(mode = "ne"))
  expect_snapshot(print(x))
})


test_that("print works for DN glycan graphs", {
  skip_on_old_win()
  x <- new_dn_glycan_graph(n_glycan_core(mode = "dn"))
  expect_snapshot(print(x))
})


test_that("print works for NE glycan graphs with verbose = FALSE", {
  x <- new_ne_glycan_graph(n_glycan_core(mode = "ne"))
  expect_snapshot(print(x, verbose = FALSE))
})


test_that("print works for DN glycan graphs with verbose = FALSE", {
  x <- new_dn_glycan_graph(n_glycan_core(mode = "dn"))
  expect_snapshot(print(x, verbose = FALSE))
})


test_that("print works for NE glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_ne_glycan_graph(n_glycan_core(mode = "ne", linkage = FALSE))
  expect_snapshot(print(x))
})


test_that("print works for DN glycan graphs without linkages", {
  skip_on_old_win()
  x <- new_dn_glycan_graph(n_glycan_core(mode = "dn", linkage = FALSE))
  expect_snapshot(print(x))
})


test_that("print works for one-mono NE glycan graphs", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_ne_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for one-mono DN glycan graphs", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  igraph::V(graph)$linkage <- NA_character_
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_dn_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for NE graphs with substituent", {
  skip_on_old_win()
  glycan <- o_glycan_core_1(mode = "ne")
  igraph::V(glycan)$sub <- "6S"
  expect_snapshot(print(glycan))
})


test_that("print works for DN graphs with substituent", {
  skip_on_old_win()
  glycan <- o_glycan_core_1(mode = "dn")
  igraph::V(glycan)$sub <- c("6S", "6S", NA)
  expect_snapshot(print(glycan))
})


test_that("print works for one-node NE graphs with substituent", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- "6S"
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_ne_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for one-node DN graphs with substituent", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- "6S"
  igraph::V(graph)$linkage <- NA_character_
  graph$anomer <- "a1"
  graph$alditol <- FALSE
  glycan <- new_dn_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for NE graphs with alditol", {
  skip_on_old_win()
  glycan <- o_glycan_core_1(mode = "ne")
  glycan$alditol <- TRUE
  expect_snapshot(print(glycan))
})


test_that("print works for DN graphs with alditol", {
  skip_on_old_win()
  glycan <- o_glycan_core_1(mode = "dn")
  glycan$alditol <- TRUE
  expect_snapshot(print(glycan))
})


test_that("print works for one-node NE graphs with alditol", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  graph$anomer <- "a1"
  graph$alditol <- TRUE
  glycan <- new_ne_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for one-node DN graphs with alditol", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- ""
  igraph::V(graph)$linkage <- NA_character_
  graph$anomer <- "a1"
  graph$alditol <- TRUE
  glycan <- new_dn_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for NE graph with alditol and substituent on root", {
  skip_on_old_win()
  glycan <- o_glycan_core_1(mode = "ne")
  igraph::V(glycan)$sub <- c("6S", "")
  glycan$alditol <- TRUE
  expect_snapshot(print(glycan))
})


test_that("print works for DN graph with alditol and substituent on root", {
  skip_on_old_win()
  glycan <- o_glycan_core_1(mode = "dn")
  igraph::V(glycan)$sub <- c("6S", "", NA)
  glycan$alditol <- TRUE
  expect_snapshot(print(glycan))
})


test_that("print works for one-node NE graph with alditol and substituent on root", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- "6S"
  graph$anomer <- "a1"
  graph$alditol <- TRUE
  glycan <- new_ne_glycan_graph(graph)
  expect_snapshot(print(glycan))
})


test_that("print works for one-node DN graph with alditol and substituent on root", {
  skip_on_old_win()
  graph <- igraph::make_graph(NULL, n = 1, directed = TRUE)
  igraph::V(graph)$type <- "mono"
  igraph::V(graph)$mono <- "N"
  igraph::V(graph)$sub <- "6S"
  igraph::V(graph)$linkage <- NA_character_
  graph$anomer <- "a1"
  graph$alditol <- TRUE
  glycan <- new_dn_glycan_graph(graph)
  expect_snapshot(print(glycan))
})
