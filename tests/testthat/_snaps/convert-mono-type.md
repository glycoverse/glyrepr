# converting glycan from simple to generic fails

    Code
      convert_mono_type(glycan, to = "generic")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `convert_mono_type()`:
      ! Cannot convert from "simple" to "generic".
      i Can only convert in this order: concrete -> generic -> simple.

# converting glycan from simple to concrete fails

    Code
      convert_mono_type(glycan, to = "concrete")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `convert_mono_type()`:
      ! Cannot convert from "simple" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting glycan from generic to concrete fails

    Code
      convert_mono_type(glycan, to = "concrete")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `convert_mono_type()`:
      ! Cannot convert from "generic" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting glycan mono types with NA produced

    Code
      convert_mono_type(modified_glycan_vec, to = "generic")
    Condition
      Error in `purrr::map()`:
      i In index: 1.
      Caused by error in `convert_mono_type()`:
      ! Some monosaccharides cannot be converted to "generic": "Pse".

# convert mono types with bad directions

    Code
      convert_mono_type(before, to)
    Condition
      Error in `convert_mono_type()`:
      ! These monosaccharides cannot be converted to "concrete": "H" and "Hex"
      i Conversion could only be done in this direction: concrete -> generic -> simple

# converting mono types with NA

    Code
      convert_mono_type(c("Pse", "Fuc"), "generic")
    Condition
      Warning in `convert_mono_type()`:
      Some monosaccharides cannot be converted to "generic": "Pse".
    Output
      [1] NA     "dHex"

# get mono types fails for unknown monos

    Code
      get_mono_type(mono)
    Condition
      Error in `get_mono_type()`:
      ! Unknown monosaccharide: "bad1", "bad2", and "bad3".

