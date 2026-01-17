# reduce_structure_level rejects higher level

    Code
      reduce_structure_level(glycan, to_level = "topological")
    Condition
      Error in `reduce_structure_level()`:
      ! Cannot reduce a structure to a higher resolution level.
      x Some structures in `x` have levels: "basic".
      i Target level: "topological" (> "basic").
      i You can use `get_structure_level()` to check the structure levels of `x`.

