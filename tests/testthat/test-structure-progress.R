test_that("progress preserves values, names, missingness and graph attributes", {
  strings <- c(a = "Glc(?1-", b = NA, c = "Glc(?1-", d = "Gal(b1-4)Glc(?1-")
  graph <- get_structure_graphs(as_glycan_structure("Glc(?1-"))
  graph$source <- "retained"
  record <- list(
    mono = "Glc",
    sub = "",
    edges = integer(),
    linkage = character(),
    anomer = "?1"
  )
  cases <- list(
    strings,
    character(),
    c(NA_character_, NA_character_),
    graph,
    list(graph, graph),
    as_glycan_structure(strings)
  )
  signature <- function(x) {
    list(
      values = as.character(x),
      graphs = lapply(attr(x, "graphs"), igraph::as_data_frame, what = "both"),
      attrs = lapply(attr(x, "graphs"), igraph::graph_attr)
    )
  }
  for (policy in c("error", "na")) {
    for (x in cases) {
      expect_equal(
        signature(as_glycan_structure(x, policy, TRUE)),
        signature(as_glycan_structure(x, policy))
      )
    }
    x <- list(a = record, missing = NULL, duplicate = record)
    expect_equal(
      signature(structure_from_arrays(x, policy, TRUE)),
      signature(structure_from_arrays(x, policy))
    )
  }
})

test_that("callbacks report native and construction stages with deduplicated totals", {
  events <- list()
  report <- function(stage, current, total) {
    events[[length(events) + 1L]] <<- list(
      stage = stage,
      current = current,
      total = total
    )
  }
  as_glycan_structure(c("Glc(?1-", NA, "Glc(?1-", "Gal(?1-"), progress = report)
  expect_equal(
    vapply(events, `[[`, character(1), "stage"),
    rep(c("Parsing unique structures", "Building structures"), each = 2)
  )
  expect_equal(vapply(events, `[[`, numeric(1), "current"), c(0, 2, 0, 2))
  expect_equal(vapply(events, `[[`, numeric(1), "total"), rep(2, 4))
  events <- list()
  record <- list(
    mono = "Glc",
    sub = "",
    edges = integer(),
    linkage = character(),
    anomer = "?1"
  )
  structure_from_arrays(list(record, NULL), progress = report)
  expect_equal(
    unique(vapply(events, `[[`, character(1), "stage")),
    c(
      "Validating records",
      "Canonicalizing records",
      "Building structures",
      "Assembling results"
    )
  )
  expect_equal(tail(events, 1)[[1]]$current, 1)
})

test_that("progress preserves failure diagnostics and recovery positions", {
  capture <- function(progress, arrays, policy) {
    warnings <- list()
    result <- withCallingHandlers(
      tryCatch(
        if (arrays) {
          structure_from_arrays(list(bad = list()), policy, progress)
        } else {
          as_glycan_structure(
            c(ok = "Glc(?1-", bad = "invalid", again = "invalid"),
            policy,
            progress
          )
        },
        error = identity
      ),
      warning = function(w) {
        warnings[[length(warnings) + 1L]] <<- w
        invokeRestart("muffleWarning")
      }
    )
    list(
      value = if (inherits(result, "error")) {
        conditionMessage(result)
      } else {
        as.character(result)
      },
      class = class(result),
      warnings = lapply(warnings, function(w) {
        list(
          message = conditionMessage(w),
          positions = w$positions,
          reasons = w$reasons
        )
      })
    )
  }
  for (arrays in c(FALSE, TRUE)) {
    for (policy in c("error", "na")) {
      expect_equal(
        capture(TRUE, arrays, policy),
        capture(FALSE, arrays, policy)
      )
    }
  }
})

test_that("callback failures are not recovered as invalid structures", {
  for (stage in c("Parsing unique structures", "Building structures")) {
    report <- function(s, current, total) {
      if (s == stage && current == total) stop("callback stopped")
    }
    error <- tryCatch(
      as_glycan_structure("Glc(?1-", "na", report),
      error = identity
    )
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "callback stopped")
  }
})

test_that("cli bars close on success, errors and interrupts without closing caller bars", {
  outer <- cli::cli_progress_bar(total = 10, .auto_close = FALSE)
  on.exit(cli::cli_progress_done(id = outer), add = TRUE)
  expect_equal(cli::cli_progress_num(), 1L)
  as_glycan_structure("Glc(?1-", progress = TRUE)
  expect_equal(cli::cli_progress_num(), 1L)
  tryCatch(as_glycan_structure("invalid", progress = TRUE), error = identity)
  expect_equal(cli::cli_progress_num(), 1L)
  local_mocked_bindings(.compact_iupac_arrays = function(...) {
    rlang::interrupt()
  })
  tryCatch(
    as_glycan_structure("Glc(?1-", progress = TRUE),
    interrupt = identity
  )
  expect_equal(cli::cli_progress_num(), 1L)
})

test_that("native callbacks cover missing and failed records without swallowing errors", {
  record <- list(
    mono = "Glc",
    sub = "",
    edges = integer(),
    linkage = character(),
    anomer = "?1"
  )
  for (stage in c(
    "Validating records",
    "Canonicalizing records",
    "Building structures",
    "Assembling results"
  )) {
    report <- function(s, current, total) {
      if (s == stage && current == total) stop("callback stopped")
    }
    error <- tryCatch(
      structure_from_arrays(list(record, NULL), "na", report),
      error = identity
    )
    expect_s3_class(error, "error")
    expect_match(conditionMessage(error), "callback stopped")
  }
  events <- list()
  report <- function(stage, current, total) {
    events[[length(events) + 1L]] <<- c(current, total)
  }
  expect_length(structure_from_arrays(list(), progress = report), 0L)
  expect_length(events, 0L)
  expect_equal(
    is.na(structure_from_arrays(list(NULL, NULL), progress = report)),
    c(TRUE, TRUE)
  )
  expect_equal(tail(events, 1)[[1]], c(2, 2))
})

test_that("graph conversion and fallback report progress", {
  stages <- character()
  report <- function(stage, current, total) stages <<- c(stages, stage)
  graph <- get_structure_graphs(as_glycan_structure("Glc(?1-"))
  as_glycan_structure(list(graph, graph), progress = report)
  expect_equal(
    unique(stages),
    c("Validating graphs", "Canonicalizing graphs", "Building structures")
  )
  stages <- character()
  local_mocked_bindings(.compact_iupac_arrays = function(x, progress = NULL) {
    lapply(x, function(x) list(status = "unsupported"))
  })
  expect_equal(
    as.character(as_glycan_structure("Glc(?1-", progress = report)),
    "Glc(?1-"
  )
  expect_equal(stages[[1]], "Parsing reference structures")
})
