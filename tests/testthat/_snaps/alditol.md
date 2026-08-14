# alditol graph attributes are validated

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `validate_glycan_graph()`:
      ! Glycan structure graph attribute alditol must be one non-missing logical value.

---

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `validate_glycan_graph()`:
      ! Glycan structure graph attribute alditol must be one non-missing logical value.

---

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `validate_glycan_graph()`:
      ! Glycan structure graph attribute alditol must be one non-missing logical value.

# alditol markers are restricted to the main reducing end

    Code
      as_glycan_structure("Gal-ol(b1-4)Glc(a1-")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `value[[3L]]()`:
      ! Could not parse IUPAC-condensed string: "Gal-ol(b1-4)Glc(a1-"
      i The alditol marker "-ol" is only allowed on the main reducing-end residue.

---

    Code
      as_glycan_structure("Glc-ol-ol(a1-")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `value[[3L]]()`:
      ! Could not parse IUPAC-condensed string: "Glc-ol-ol(a1-"
      i The alditol marker "-ol" is only allowed on the main reducing-end residue.

---

    Code
      as_glycan_structure("{Gal-ol(b1-4)|2}Glc(a1-")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! A floating glycan part cannot contain an alditol reducing end.
