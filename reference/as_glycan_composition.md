# Convert to Glycan Composition

Convert an object to a glycan composition. The resulting composition can
contain both monosaccharides and substituents.

## Usage

``` r
as_glycan_composition(x)

# S3 method for class 'glyrepr_composition'
as_glycan_composition(x)

# S3 method for class 'glyrepr_structure'
as_glycan_composition(x)

# S3 method for class 'character'
as_glycan_composition(x)

# Default S3 method
as_glycan_composition(x)
```

## Arguments

- x:

  An object to convert to a glycan composition. Can be a named integer
  vector, a list of named integer vectors, a glycan structure vector, or
  an existing glyrepr_composition object.

## Value

A glyrepr_composition object.

## Details

When converting from glycan structures, both monosaccharides and
substituents are counted. Substituents are extracted from the `sub`
attribute of each vertex in the structure. For example, a vertex with
`sub = "3Me"` contributes one "Me" substituent to the composition.

## Examples

``` r
# From a named vector
as_glycan_composition(c(Hex = 5, HexNAc = 2))
#> <glycan_composition[1]>
#> [1] Hex(5)HexNAc(2)

# From a named vector with substituents
as_glycan_composition(c(Glc = 2, Gal = 1, Me = 1, S = 1))
#> <glycan_composition[1]>
#> [1] Glc(2)Gal(1)Me(1)S(1)

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
as_glycan_composition(c("H5N2", "H5N4S1F1", "H5N4A1G1"))
#> <glycan_composition[3]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(5)HexNAc(4)dHex(1)NeuAc(1)
#> [3] Hex(5)HexNAc(4)NeuAc(1)NeuGc(1)

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

# From structures with substituents
# (This will count both monosaccharides and any substituents present)
```
