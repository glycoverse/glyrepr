# converting glycan from simple to generic fails

    Code
      convert_glycan_mono_type(glycan, to = "generic")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! Cannot convert from "simple" to "generic".
      i Can only convert in this order: concrete -> generic -> simple.

# converting glycan from simple to concrete fails

    Code
      convert_glycan_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! Cannot convert from "simple" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting glycan from generic to concrete fails

    Code
      convert_glycan_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! Cannot convert from "generic" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting glycan from generic to generic fails

    Code
      convert_glycan_mono_type(glycan, to = "generic")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! It is already "generic".

# convert mono types fails for monos already in simple form

    Code
      convert_mono_type(before, to)
    Condition
      Error in `convert_mono_type()`:
      ! These monosaccharides are already "simple": "H", "H", and "N"

# convert mono types with bad directions

    Code
      convert_mono_type(before, to)
    Condition
      Error in `convert_mono_type()`:
      ! These monosaccharides cannot be converted to "concrete": "H" and "Hex"
      i Conversion could only be done in this direction: concrete -> generic -> simple

# deciding mono types fails for multiple monos

    Code
      decide_mono_type(mono)
    Condition
      Error in `decide_mono_type()`:
      ! Unknown monosaccharide: "bad1", "bad2", and "bad3".

