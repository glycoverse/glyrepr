# malformed floating-part IUPAC is rejected

    Code
      cat(paste0(invalid, "\n", purrr::map_chr(errors, conditionMessage), collapse = "\n\n"))
    Output
      {Neu5Ac(a2-3)}
      i In index: 1.
      Caused by error in `split_floating_iupac()`:
      ! A floating glycan structure must have a main glycan.
      
      {}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! A floating part cannot be empty.
      
      {Neu5Ac}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! A floating part must end with its linkage to the main glycan.
      
      {Neu5Ac(a2-3)|}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parents must be comma-separated positive node indices.
      
      {Neu5Ac(a2-3)|0}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parents must be comma-separated positive node indices.
      
      {Neu5Ac(a2-3)|1,1}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parent indices must be unique.
      
      {Neu5Ac(a2-3)|3}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part parent index is outside the main glycan.
      i The main glycan has 2 nodes.
      
      {Neu5Ac(a2-3)|1|2}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! A floating part can contain at most one parent-index separator.
      
      {{Neu5Ac(a2-3)}}Gal(b1-3)GalNAc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac()`:
      ! Floating parts cannot be nested.
      
      Gal(b1-3)GalNAc(a1-{Neu5Ac(a2-3)}
      i In index: 1.
      Caused by error in `value[[3L]]()`:
      ! Could not parse IUPAC-condensed string: "Gal(b1-3)GalNAc(a1-{Neu5Ac(a2-3)}"
      i Invalid characters or format in IUPAC-condensed string

# invalid floating-part graph annotations are rejected

    Code
      cat(paste0(names(errors), "\n", purrr::map_chr(errors, conditionMessage),
      collapse = "\n\n"))
    Output
      unannotated
      i In index: 1.
      Caused by error in `validate_floating_graph_shape()`:
      ! Glycan structure must be an out tree.
      
      malformed_attr
      i In index: 1.
      Caused by error in `normalize_floating_parts()`:
      ! Graph attribute floating_parts must be a list.
      
      nonroot
      i In index: 1.
      Caused by error in `validate_floating_graph_shape()`:
      ! Each floating part root must be the root of its graph component.
      
      floating_parent
      i In index: 1.
      Caused by error in `validate_floating_graph_shape()`:
      ! Floating part parent indices must refer to vertices in the main tree.
      
      duplicate_parents
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part parent indices must be unique.
      
      duplicate_component
      i In index: 1.
      Caused by error in `floating_graph_info()`:
      ! Each floating part must identify a different graph component.
      
      undeclared_component
      i In index: 1.
      Caused by error in `floating_graph_info()`:
      ! A floating glycan structure must contain exactly one main tree.

# floating parent indices reject integer overflow

    Code
      as_glycan_structure("{Neu5Ac(a2-3)|2147483648}Gal(a1-")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parent indices exceed the supported integer range.

---

    Code
      glycan_structure(graph)
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part parents must contain valid vertex indices.

