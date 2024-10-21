# converting from simple to generic fails

    Code
      convert_glycan_mono_type(glycan, to = "generic")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! Cannot convert from "simple" to "generic".
      i Can only convert in this order: concrete -> generic -> simple.

# converting from simple to concrete fails

    Code
      convert_glycan_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! Cannot convert from "simple" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting from generic to concrete fails

    Code
      convert_glycan_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! Cannot convert from "generic" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting from generic to generic fails

    Code
      convert_glycan_mono_type(glycan, to = "generic")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! It is already "generic".

# converting from simple to simple fails

    Code
      convert_glycan_mono_type(glycan, to = "simple")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! It is already "simple".

# converting from concrete to concrete fails

    Code
      convert_glycan_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_glycan_mono_type()`:
      ! It is already "concrete".

