# malformed array records fail before graph construction

    Code
      structure_from_arrays(list(list(mono = "Gal")))
    Condition
      Error in `structure_from_arrays()`:
      ! Invalid structure at position 1.
      Caused by error in `.validate_structure_arrays()`:
      ! Structure records require "mono", "sub", "edges", "linkage", and "anomer".

---

    Code
      structure_from_arrays(list(a))
    Condition
      Error in `structure_from_arrays()`:
      ! Invalid structure at position 1.
      Caused by error in `.check_array_indices()`:
      ! Assertion on 'x' failed: Element 2 is not <= 1.

---

    Code
      result <- structure_from_arrays(list(bad = a, absent = NULL), on_failure = "na")
    Condition
      Warning:
      1 structure failed validation and was replaced with `NA`.
      x Position 1 (`bad`): Assertion on 'x' failed: Element 2 is not <= 1.

# native topology validation rejects cycles and incomplete components

    Code
      result <- structure_from_arrays(bad, on_failure = "na")
    Condition
      Warning:
      3 structures failed validation and were replaced with `NA`.
      x Position 1: Glycan structure must be an out tree., Position 2: Glycan structure must be an out tree., and Position 3: Floating part node metadata does not match its graph component. x Floating part 1 must contain exactly node 2.

