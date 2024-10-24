patrick::with_parameters_test_that("`valid_linkages()` works for good linkages", {
  expect_true(valid_linkages(char))
}, char = c("b1-4", "a1-3"))


patrick::with_parameters_test_that("`valid_linkages()` works for known linkages", {
  expect_true(valid_linkages(char))
}, char = c("b1-?", "a?-3", "?1-4", "??-?"))


patrick::with_parameters_test_that("`valid_linkages()` fails for bad linkages", {
  expect_false(valid_linkages(char))
}, char = c("1-4", "c1-4", "b1", "abc", ""))
