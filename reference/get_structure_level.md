# Get the Structure Resolution Levels

Glycan structures can have four possible levels of resolution:

- "intact": All monosaccharides are concrete (e.g. "Man", "GlcNAc"), and
  all linkages are fully determined (e.g. "a2-3", "b1-4").

- "partial": All monosaccharides are concrete (e.g. "Man", "GlcNAc"),
  but some linkage information is missing (e.g. "a2-?").

- "topological": All monosaccharides are concrete (e.g. "Man",
  "GlcNAc"), but the linkage information is completely unknown ("??-?").

- "basic": All monosaccharides are generic (e.g. "Hex", "HexNAc"), and
  the linkage information is completely unknown ("??-?").

Note that in theory you can have a glycan with generic monosaccharides
with all linkages determined. For example, "Hex(b1-3)HexNAc(a1-" is a
valid glycan structure. But in reality, this is almost impossible,
because linkage information is far more difficult to acquire than
monosaccharide information. This kind of glycan structure is also
assigned to "basic" level.

## Usage

``` r
get_structure_level(x)
```

## Arguments

- x:

  A
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  vector.

## Value

A character vector of the same length as `x`, containing the structure
level for each element.

## See also

[`has_linkages()`](https://glycoverse.github.io/glyrepr/reference/has_linkages.md),
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/reference/get_mono_type.md)

## Examples

``` r
structures <- as_glycan_structure(c(
  "Gal(b1-3)GalNAc(a1-",
  "Gal(b1-?)GalNAc(a1-",
  "Gal(??-?)GalNAc(??-",
  "Hex(??-?)HexNAc(??-",
  "Hex(b1-3)HexNAc(a1-"
))
get_structure_level(structures)
#> [1] "intact"      "partial"     "topological" "basic"       "basic"      
```
