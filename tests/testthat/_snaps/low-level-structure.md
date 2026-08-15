# empty floating metadata preserves the validation contract

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `normalize_floating_parts()`:
      ! Graph attribute floating_parts must be a list.

# validate_glycan_graph reports invalid input

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `validate_glycan_graph()`:
      ! Glycan structure must be directed.

# graph_to_iupac generates one string from one graph

    Code
      graph_to_iupac(list(graph))
    Condition
      Error in `graph_to_iupac()`:
      ! Assertion on 'graph' failed: Must inherit from class 'igraph', but has class 'list'.

# unannotated forests remain invalid

    Code
      validate_glycan_graph(graph)
    Condition
      Error in `validate_floating_graph_shape()`:
      ! Glycan structure must be an out tree.

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

