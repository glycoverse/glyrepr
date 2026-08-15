# Get the Structure Resolution Levels

Glycan structures can have three possible levels of resolution,
determined only by linkage and anomer information:

- "intact": No linkage or anomer is unknown or ambiguous.

- "partial": At least one linkage or anomer is unknown or ambiguous, and
  at least one linkage or anomer contains known information.

- "topological": All linkages and anomers are completely unknown
  ("??-?"/"??").

Floating metadata does not by itself affect the structure level. A
floating part's attachment linkage is evaluated like an ordinary
linkage, while candidate-parent ambiguity for floating parts and
substituents is independent of linkage resolution. Consequently, a
glycan with floating metadata is regarded as "intact" when every
graph-edge and floating-part attachment linkage, as well as the
reducing-end anomer, is fully specified, even though a parent residue
remains unlocalized.

## Usage

``` r
get_structure_level(x)
```

## Arguments

- x:

  A
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
  vector or a glycan `igraph`.

## Value

For vector input, a character vector containing one structure level per
element. Missing elements return `NA_character_`, and an empty input
returns [`character()`](https://rdrr.io/r/base/character.html). For
graph input, a character scalar.

## See also

[`has_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/has_linkages.md),
[`has_floating_parts()`](https://glycoverse.github.io/glyrepr/dev/reference/has_floating_parts.md),
[`has_floating_substituents()`](https://glycoverse.github.io/glyrepr/dev/reference/has_floating_substituents.md),
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md)

## Examples

``` r
glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
get_structure_level(glycan)
#> [1] "intact"

floating <- as_glycan_structure(
  "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
)
get_structure_level(floating)
#> [1] "intact"
```
