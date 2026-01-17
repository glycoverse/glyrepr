# Create a Glycan Composition

Create a glycan composition from a list of named integer vectors.
Compositions can contain both monosaccharides and substituents.

## Usage

``` r
glycan_composition(...)

is_glycan_composition(x)
```

## Arguments

- ...:

  Named integer vectors. Names are monosaccharides or substituents,
  values are numbers of residues. Monosaccharides and substituents can
  be mixed in the same composition.

- x:

  A list of named integer vectors.

## Value

A glyrepr_composition object.

## Details

Compositions can contain:

- Monosaccharides: either generic (e.g., "Hex", "HexNAc") or concrete
  (e.g., "Glc", "Gal"). All monosaccharides in a composition vector must
  be of the same type.

- Substituents: e.g., "Me", "Ac", "S". These can be mixed with either
  generic or concrete monosaccharides.

Components are automatically sorted with monosaccharides first
(according to their order in the monosaccharides table), followed by
substituents (according to their order in
[`available_substituents()`](https://glycoverse.github.io/glyrepr/reference/available_substituents.md)).

## See also

[`available_monosaccharides()`](https://glycoverse.github.io/glyrepr/reference/available_monosaccharides.md),
[`available_substituents()`](https://glycoverse.github.io/glyrepr/reference/available_substituents.md)

## Examples

``` r
# A vector with one composition (generic monosaccharides)
glycan_composition(c(Hex = 5, HexNAc = 2))
#> <glycan_composition[1]>
#> [1] Hex(5)HexNAc(2)
# A vector with multiple compositions
glycan_composition(c(Hex = 5, HexNAc = 2), c(Hex = 5, HexNAc = 4, dHex = 2))
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(5)HexNAc(4)dHex(2)
# Residues are reordered automatically
glycan_composition(c(HexNAc = 1, Hex = 2))
#> <glycan_composition[1]>
#> [1] Hex(2)HexNAc(1)
# An example for generic monosaccharides
glycan_composition(c(Hex = 2, HexNAc = 1))
#> <glycan_composition[1]>
#> [1] Hex(2)HexNAc(1)
# An example for concrete monosaccharides
glycan_composition(c(Glc = 2, Gal = 1))
#> <glycan_composition[1]>
#> [1] Glc(2)Gal(1)
# Compositions with substituents
glycan_composition(c(Glc = 1, S = 1))
#> <glycan_composition[1]>
#> [1] Glc(1)S(1)
glycan_composition(c(Hex = 3, HexNAc = 2, Me = 1, Ac = 1))
#> <glycan_composition[1]>
#> [1] Hex(3)HexNAc(2)Me(1)Ac(1)
# Substituents are sorted after monosaccharides
glycan_composition(c(S = 1, Gal = 1, Ac = 1, Glc = 1))
#> <glycan_composition[1]>
#> [1] Glc(1)Gal(1)Ac(1)S(1)
```
