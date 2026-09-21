# batch graph recovery distinguishes missing and invalid inputs

    Code
      result <- canonicalize_glycan_graphs(list(absent = NULL, bad = 1), on_failure = "na")
    Condition
      Warning:
      1 structure failed validation and was replaced with `NA`.
      x Position 2 (`bad`): Assertion on 'graph' failed: Must inherit from class 'igraph', but has class 'numeric'.

---

    Code
      canonicalize_glycan_graphs(list(NULL, 1))
    Condition
      Error in `canonicalize_glycan_graphs()`:
      ! Invalid structure at position 2.
      Caused by error in `doTryCatch()`:
      ! Assertion on 'graph' failed: Must inherit from class 'igraph', but has class 'numeric'.

