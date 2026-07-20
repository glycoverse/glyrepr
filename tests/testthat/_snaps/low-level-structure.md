# validate_glycan_graph reports invalid input

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `validate_glycan_graph()`:
      ! Glycan structure must be directed.

# validate_glycan_graph_vector rejects mixed monosaccharide types

    Code
      validate_glycan_graph_vector(list(concrete, generic))
    Condition
      Error in `validate_glycan_graph_vector()`:
      ! All structures must have the same monosaccharide type.
      x Found 1 concrete and 1 generic structure(s) in the same vector.
      i Use `convert_to_generic()` to convert concrete structures to generic type.

# graph_to_iupac generates one string from one graph

    Code
      graph_to_iupac(list(graph))
    Condition
      Error in `graph_to_iupac()`:
      ! Assertion on 'graph' failed: Must inherit from class 'igraph', but has class 'list'.

# new_glycan_structure checks graph lookup integrity

    Code
      new_glycan_structure(iupac, list(graph))
    Condition
      Error in `new_glycan_structure()`:
      ! `graphs` must have unique, non-missing IUPAC-condensed names.

---

    Code
      new_glycan_structure("missing-key", stats::setNames(list(graph), iupac))
    Condition
      Error in `new_glycan_structure()`:
      ! `graphs` does not contain every structure in `iupac`.
      x Missing graph: "missing-key".

