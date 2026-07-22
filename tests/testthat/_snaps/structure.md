# as_glycan_structure can replace invalid graphs with NA

    Code
      result <- as_glycan_structure(graphs, on_failure = "na")
    Condition
      Warning:
      1 structure failed validation and was replaced with `NA`.
      x Position 3 (`invalid`): Unknown monosaccharide: NotAMonosaccharide

# as_glycan_structure can replace invalid character input with NA

    Code
      result <- as_glycan_structure(iupacs, on_failure = "na")
    Condition
      Warning:
      1 structure failed validation and was replaced with `NA`.
      x Position 2 (`invalid`): Could not parse IUPAC-condensed string: "not-a-structure" i Can't extract anomer information. i Anomer information is required for the reducing-end monosaccharide. i For example, use 'Man(a1-' instead of 'Man'.

# as_glycan_structure keeps strict failures as the default

    Code
      as_glycan_structure(list(valid, invalid))
    Condition
      Error in `purrr::map()`:
      i In index: 2.
      Caused by error in `validate_glycan_graph()`:
      ! Unknown monosaccharide: NotAMonosaccharide

# as_glycan_structure keeps vector-level failures strict

    Code
      as_glycan_structure(iupacs, on_failure = "na")
    Condition
      Error in `validate_glycan_graph_vector()`:
      ! All structures must have the same monosaccharide type.
      x Found 1 concrete and 1 generic structure(s) in the same vector.
      i Use `convert_to_generic()` to convert concrete structures to generic type.

# as_glycan_structure validates on_failure

    Code
      as_glycan_structure("Glc(?1-", on_failure = "skip")
    Condition
      Error in `as_glycan_structure()`:
      ! `on_failure` must be one of "error" or "na", not "skip".

# format.glyrepr_structure includes names with tab separation

    Code
      format(glycans)
    Output
      [1] "A\tGal(b1-3)GalNAc(a1-                                "
      [2] "B\tMan(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

# format.glyrepr_structure without names works correctly

    Code
      format(glycans)
    Output
      [1] "Gal(b1-3)GalNAc(a1-                                "
      [2] "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"

# truncation works in tibble

    Code
      print(tibble, width = 30)
    Output
      # A tibble: 3 x 2
        struc                      a
        <struct>               <dbl>
      1 Man(a1-3)[Man(a1-6)]M~     1
      2 Man(a1-3)[Man(a1-6)]M~     1
      3 Man(a1-3)[Man(a1-6)]M~     1

# print.glyrepr_structure supports n

    <glycan_structure[11]>
    [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [2] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [3] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [4] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [5] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [6] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [7] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [8] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [9] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [10] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    ... (1 more not shown)
    # Unique structures: 1

---

    <glycan_structure[11]>
    [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [2] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [3] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [4] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [5] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [6] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [7] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [8] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [9] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [10] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    [11] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
    # Unique structures: 1

