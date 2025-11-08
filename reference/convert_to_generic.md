# Convert Monosaccharides to Generic Type

This function converts monosaccharide types of monosaccharide
characters, glycan compositions, or glycan structures from concrete to
generic type. This is a simplified version that only supports conversion
from "concrete" to "generic" monosaccharides.

## Usage

``` r
convert_to_generic(x)

# S3 method for class 'character'
convert_to_generic(x)

# S3 method for class 'glyrepr_structure'
convert_to_generic(x)

# S3 method for class 'glyrepr_composition'
convert_to_generic(x)
```

## Arguments

- x:

  Either of these objects:

  - A character of monosaccharide;

  - A glycan composition vector ("glyrepr_composition" object);

  - A glycan structure vector ("glyrepr_structure" object).

## Value

A new object of the same class as `x` with monosaccharides converted to
generic type.

## Two types of monosaccharides

There are two types of monosaccharides:

- concrete: e.g. "Gal", "GlcNAc", "Glc", "Fuc", etc.

- generic: e.g. "Hex", "HexNAc", "HexA", "HexN", etc.

For the full list of monosaccharides, use
[`available_monosaccharides()`](https://glycoverse.github.io/glyrepr/reference/available_monosaccharides.md).

## Examples

``` r
# Convert character vectors
convert_to_generic(c("Gal", "GlcNAc"))
#> [1] "Hex"    "HexNAc"

# Convert glycan compositions
comps <- glycan_composition(
  c(Gal = 5, GlcNAc = 2),
  c(Glc = 5, GalNAc = 4, Fuc = 1)
)
convert_to_generic(comps)
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(5)HexNAc(4)dHex(1)

# Convert glycan structures
strucs <- glycan_structure(
  n_glycan_core(),
  o_glycan_core_1()
)
convert_to_generic(strucs)
#> <glycan_structure[2]>
#> [1] Hex(a1-3)[Hex(a1-6)]Hex(b1-4)HexNAc(b1-4)HexNAc(b1-
#> [2] Hex(b1-3)HexNAc(a1-
#> # Unique structures: 2
```
