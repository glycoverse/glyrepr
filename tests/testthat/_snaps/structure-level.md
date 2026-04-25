# get_structure_level works for a basic glycan vector with linkages

    Code
      res <- get_structure_level(glycans)
    Condition
      Warning:
      Generic glycan structures with linkage annotations are treated as "basic".
      i Linkage information is ignored when residues are generic.

# reduce_structure_level rejects higher level

    Code
      reduce_structure_level(glycan, to_level = "topological")
    Condition
      Error in `reduce_structure_level()`:
      ! Cannot reduce a structure to a higher resolution level.
      x Structure level of `x`: "basic".
      i Target level: "topological" (> "basic").
      i You can use `get_structure_level()` to check the structure level of `x`.

