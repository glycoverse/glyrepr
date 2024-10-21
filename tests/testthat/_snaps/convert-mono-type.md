# converting from simple to generic fails

    Code
      convert_mono_type(glycan, to = "generic")
    Condition
      Error in `convert_mono_type()`:
      ! Cannot convert from "simple" to "generic".
      i Can only convert in this order: concrete -> generic -> simple.

# converting from simple to concrete fails

    Code
      convert_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_mono_type()`:
      ! Cannot convert from "simple" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

# converting from generic to concrete fails

    Code
      convert_mono_type(glycan, to = "concrete")
    Condition
      Error in `convert_mono_type()`:
      ! Cannot convert from "generic" to "concrete".
      i Can only convert in this order: concrete -> generic -> simple.

