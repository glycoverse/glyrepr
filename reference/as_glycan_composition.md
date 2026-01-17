# Convert to Glycan Composition

Converts an object to a glycan composition using vctrs casting
framework. This function provides a convenient way to convert various
input types to
[`glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.md).

## Usage

``` r
as_glycan_composition(x)
```

## Arguments

- x:

  An object to convert to a glycan composition. Supported inputs
  include:

  - Named integer vectors or lists of named integer vectors

  - Character vectors with composition strings (e.g., "Hex(5)HexNAc(2)")

  - `glyrepr_structure` objects (counts both monosaccharides and
    substituents)

  - Existing `glyrepr_composition` objects (returned as-is)

## Value

A `glyrepr_composition` object.

## Details

This function uses the vctrs casting framework for type conversion. When
converting from glycan structures, both monosaccharides and substituents
are counted. Substituents are extracted from the `sub` attribute of each
vertex in the structure. For example, a vertex with `sub = "3Me"`
contributes one "Me" substituent to the composition.

## Examples

``` r
# From a single named vector
as_glycan_composition(c(Hex = 5, HexNAc = 2))
#> <glycan_composition[1]>
#> [1] Hex(5)HexNAc(2)

# From a list of named vectors
as_glycan_composition(list(c(Hex = 5, HexNAc = 2), c(Hex = 3, HexNAc = 1)))
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(3)HexNAc(1)

# From a character vector of Byonic composition strings
as_glycan_composition(c("Hex(5)HexNAc(2)", "Hex(3)HexNAc(1)"))
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(3)HexNAc(1)

# From a character vector of simple composition strings
as_glycan_composition(c("H5N2", "H5N4S1F1"))
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(5)HexNAc(4)dHex(1)NeuAc(1)

# From an existing composition (returns as-is)
comp <- glycan_composition(c(Hex = 5, HexNAc = 2))
as_glycan_composition(comp)
#> <glycan_composition[1]>
#> [1] Hex(5)HexNAc(2)

# From a glycan structure vector
strucs <- c(n_glycan_core(), o_glycan_core_1())
as_glycan_composition(strucs)
#> <glycan_composition[2]>
#> [1] Man(3)GlcNAc(2)
#> [2] Gal(1)GalNAc(1)
```
