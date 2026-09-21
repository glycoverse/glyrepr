# character recovery preserves repeated failures and canonical duplicates

    Code
      result <- as_glycan_structure(x, on_failure = "na")
    Condition
      Warning:
      2 structures failed validation and were replaced with `NA`.
      x Position 2 (`bad`): Could not parse IUPAC-condensed string: "not-a-structure" i Invalid characters or format in IUPAC-condensed string and Position 4 (`bad_again`): Could not parse IUPAC-condensed string: "not-a-structure" i Invalid characters or format in IUPAC-condensed string

# strict character errors retain parse-before-validation precedence

    Code
      as_glycan_structure(c("Gal(b1-4)[Man(a1-4)]Glc", "not-a-structure"))
    Condition
      Error in `purrr::map()`:
      i In index: 2.
      Caused by error in `value[[3L]]()`:
      ! Could not parse IUPAC-condensed string: "not-a-structure"
      i Invalid characters or format in IUPAC-condensed string

