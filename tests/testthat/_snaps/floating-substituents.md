# floating substituent assignments must be conflict-free

    Code
      cat(paste0(invalid, "\n", purrr::map_chr(errors, conditionMessage), collapse = "\n\n"))
    Output
      {6S|1}Gal6Me(a1-
      i In index: 1.
      Caused by error in `floating_substituent_domain()`:
      ! Floating substituent has impossible explicit parent metadata.
      x Floating substituent 1 cannot use explicit parent node 1 because every declared carbon position is already occupied.
      
      {6S|1,2}Gal6Me(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `floating_substituent_domain()`:
      ! Floating substituent has impossible explicit parent metadata.
      x Floating substituent 1 cannot use explicit parent node 1 because every declared carbon position is already occupied.
      
      {6S}{6Ac}Gal(a1-
      i In index: 1.
      Caused by error in `validate_floating_substituent_slots()`:
      ! Floating substituents cannot be localized simultaneously.
      x No conflict-free assignment exists for the declared parent residues and carbon positions.
      
      {4/6S|1}{4Ac|1}{6Me|1}Gal(a1-
      i In index: 1.
      Caused by error in `validate_floating_substituent_slots()`:
      ! Floating substituents cannot be localized simultaneously.
      x No conflict-free assignment exists for the declared parent residues and carbon positions.

# malformed floating substituent IUPAC is rejected

    Code
      cat(paste0(invalid, "\n", purrr::map_chr(errors, conditionMessage), collapse = "\n\n"))
    Output
      {S}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! A floating part must end with its linkage to the main glycan.
      
      {0S}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! A floating part must end with its linkage to the main glycan.
      
      {6X}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! A floating part must end with its linkage to the main glycan.
      
      {6S|}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parents must be comma-separated positive node indices.
      
      {6S|0}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parents must be comma-separated positive node indices.
      
      {6S|1,1}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! Floating part parent indices must be unique.
      
      {6S|3}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating part parent index is outside the main glycan.
      i The main glycan has 2 nodes.
      
      {6S|1|2}Gal(a1-3)Glc(a1-
      i In index: 1.
      Caused by error in `split_floating_iupac_part()`:
      ! A floating part can contain at most one parent-index separator.

# invalid floating substituent graph metadata is rejected

    Code
      cat(paste0(names(errors), "\n", purrr::map_chr(errors, conditionMessage),
      collapse = "\n\n"))
    Output
      malformed_attr
      i In index: 1.
      Caused by error in `normalize_floating_substituents()`:
      ! Graph attribute floating_substituents must be a list.
      
      missing_field
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating substituent metadata is incomplete.
      x Missing field: parents.
      
      unknown_substituent
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Each floating substituent must contain one valid token in substituent.
      
      noncanonical_substituent
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating substituent field substituent must use canonical position order.
      
      duplicate_parents
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating substituent parent indices must be unique.
      
      invalid_parent
      i In index: 1.
      Caused by error in `purrr::map()`:
      i In index: 1.
      Caused by error in `.f()`:
      ! Floating substituent parent indices must refer to valid graph vertices.
      
      floating_component_parent
      i In index: 1.
      Caused by error in `validate_floating_substituent_parents()`:
      ! Floating substituent parent indices must refer to vertices in the main tree.

