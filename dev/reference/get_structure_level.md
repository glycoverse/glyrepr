# Get the Structure Resolution Levels

Glycan structures can have four possible levels of resolution:

- "intact": All monosaccharides are concrete (e.g. "Man", "GlcNAc"), and
  no linkage or anomer is unknown or ambiguous.

- "partial": All monosaccharides are concrete (e.g. "Man", "GlcNAc"), at
  least one linkage or anomer is unknown or ambiguous, and at least one
  linkage or anomer has a non-"?" annotation.

- "topological": All monosaccharides are concrete (e.g. "Man",
  "GlcNAc"), and all linkages and anomers are completely unknown
  ("??-?"/"??").

- "basic": All monosaccharides are generic (e.g. "Hex", "HexNAc").

Floating parts do not by themselves affect the structure level. Their
attachment linkage is evaluated like an ordinary linkage, while
candidate-parent ambiguity is independent of linkage resolution.
Consequently, a glycan with floating parts is regarded as "intact" when
its residues are concrete and every graph-edge and floating-attachment
linkage, as well as the reducing-end anomer, is fully specified, even
though a floating part's parent remains unlocalized.

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
[`has_floating_parts()`](https://glycoverse.github.io/glyrepr/dev/reference/has_floating_parts.md),
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md)

## Examples

``` r
glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
get_structure_level(glycan)
#> [1] "intact"

floating <- as_glycan_structure(
  "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
)
get_structure_level(floating)
#> [1] "intact"
```
