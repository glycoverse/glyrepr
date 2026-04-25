# Get the Structure Resolution Levels

Glycan structures can have four possible levels of resolution:

- "intact": All monosaccharides are concrete (e.g. "Man", "GlcNAc"), and
  no linkage or anomer contains "?".

- "partial": All monosaccharides are concrete (e.g. "Man", "GlcNAc"), at
  least one linkage or anomer contains "?", and at least one linkage or
  anomer has a non-"?" annotation.

- "topological": All monosaccharides are concrete (e.g. "Man",
  "GlcNAc"), and all linkages and anomers are completely unknown
  ("??-?"/"??").

- "basic": All monosaccharides are generic (e.g. "Hex", "HexNAc").

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
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
  vector.

## Value

A character scalar containing the structure level for `x`. If `x` is
empty or all structures in `x` are NA, returns NA_character\_.

## See also

[`has_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/has_linkages.md),
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md)

## Examples

``` r
glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
get_structure_level(glycan)
#> [1] "intact"
```
