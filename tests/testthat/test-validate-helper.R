patrick::with_parameters_test_that(
  "`valid_linkages()` works for good linkages",
  {
    expect_true(valid_linkages(char))
  },
  char = c("b1-4", "a1-3")
)


patrick::with_parameters_test_that(
  "`valid_linkages()` works for unknown linkages",
  {
    expect_true(valid_linkages(char))
  },
  char = c("b1-?", "a?-3", "?1-4", "??-?", "a2-3/6")
)


patrick::with_parameters_test_that(
  "`valid_linkages()` fails for bad linkages",
  {
    expect_false(valid_linkages(char))
  },
  char = c("1-4", "c1-4", "a3-1", "a1-10", "b1", "abc", "")
)


test_that("linkage pattern helper matches valid linkage grammar", {
  pattern <- linkage_pattern()

  expect_true(all(stringr::str_detect(
    c("a1-2", "?1-4", "a?-?", "a2-3/6"),
    pattern
  )))
  expect_false(any(stringr::str_detect(
    c("1-4", "c1-4", "a3-1", "a1-10"),
    pattern
  )))
})
