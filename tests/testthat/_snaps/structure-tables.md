# structure_from_tibbles validates alditol status

    Code
      structure_from_tibbles(nodes, edges, "a1", alditols = logical())
    Condition
      Error in `validate_structure_alditols()`:
      ! `alditols` must have length one or the same length as `anomers`.

---

    Code
      structure_from_tibbles(nodes, edges, "a1", alditols = NA)
    Condition
      Error in `validate_structure_alditols()`:
      ! `alditols` cannot be missing for a non-missing glycan.

---

    Code
      structure_from_tibbles(nodes, edges, "a1", alditols = "FALSE")
    Condition
      Error in `validate_structure_alditols()`:
      ! Assertion on 'alditols' failed: Must be of type 'logical', not 'character'.
