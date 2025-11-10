# Get Monosaccharide Types

This function determines the type of monosaccharides in character
vectors, glycan compositions, or glycan structures. Supported types:
"concrete" and "generic" (see details below).

## Usage

``` r
get_mono_type(x)

# S3 method for class 'character'
get_mono_type(x)

# S3 method for class 'glyrepr_structure'
get_mono_type(x)

# S3 method for class 'glyrepr_composition'
get_mono_type(x)
```

## Arguments

- x:

  Either of these objects:

  - A character vector of monosaccharide names;

  - A glycan composition vector ("glyrepr_composition" object);

  - A glycan structure vector ("glyrepr_structure" object).

## Value

A character vector specifying the monosaccharide type(s). For structures
and compositions, returns the type for each element.

## Two types of monosaccharides

There are two types of monosaccharides:

- concrete: e.g. "Gal", "GlcNAc", "Glc", "Fuc", etc.

- generic: e.g. "Hex", "HexNAc", "HexA", "HexN", etc.

For the full list of monosaccharides, use
[`available_monosaccharides()`](https://glycoverse.github.io/glyrepr/reference/available_monosaccharides.md).

## Special monosaccharides

Some monosaccharides have the same name for both generic and concrete
types. For example, "Mur" is both a generic and concrete monosaccharide.
`get_mono_type()` will return NA for these monosaccharides.

If a
[`glycan_composition()`](https://glycoverse.github.io/glyrepr/reference/glycan_composition.md)
or
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
contains these special monosaccharides, the logic will be as follows:

- If at least one residue is not a special monosaccharide, the type will
  be determined as the type of the residues.

- If all residues are special monosaccharides, the type will be
  determined as "concrete".

## See also

[`convert_to_generic()`](https://glycoverse.github.io/glyrepr/reference/convert_to_generic.md)

## Examples

``` r
# Character vector
get_mono_type(c("Gal", "Hex"))
#> [1] "concrete" "generic" 

# Glycan structures
get_mono_type(n_glycan_core(mono_type = "concrete"))
#> [1] "concrete"
get_mono_type(n_glycan_core(mono_type = "generic"))
#> [1] "generic"

# Glycan compositions
comp <- glycan_composition(c(Glc = 2, GalNAc = 1))
get_mono_type(comp)
#> [1] "concrete"

# Special cases
comps <- glycan_composition(
  c(Neu = 1),
  c(Neu = 1, Glc = 1),
  c(Mur = 1, Hex = 1),
)
get_mono_type(comps)
#> [1] "concrete" "concrete" "generic" 
```
