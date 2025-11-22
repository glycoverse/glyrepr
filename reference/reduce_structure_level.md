# Reduce a Glycan Structure to a Lower Resolution Level

This function reduces a glycan structure from a higher resolution level
to a lower resolution level (see
[`get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.md)
for four possible levels of resolution). For example, it can reduce an
"intact" structure to a "topological" structure, or a "partial"
structure to a "basic" structure. One exception is that you can never
reduce an "intact" structure to "partial" level, because the "partial"
level is not deterministic.

## Usage

``` r
reduce_structure_level(x, to_level)
```

## Arguments

- x:

  A
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  vector.

- to_level:

  The resolution level to reduce to. Can be "basic" or "topological".
  Must be a lower resolution level than any structure in `x` ("intact"
  \> "partial" \> "topological" \> "basic"). If `to_level` is the same
  as some structure in `x`, the result will be the same as the input.
  You can use
  [`get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.md)
  to check the structure levels of `x`.

## Value

A
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
vector reduced to the given resolution level.

## Details

The logic is as follows:

- If `to_level` is "topological", this function calls
  [`remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.md)
  to remove all linkages.

- If `to_level` is "basic", this function calls
  [`remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.md)
  to remove all linkages, and
  [`convert_to_generic()`](https://glycoverse.github.io/glyrepr/reference/convert_to_generic.md)
  to convert all monosaccharides to generic.

## See also

[`get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.md)

## Examples

``` r
glycan <- as_glycan_structure("Gal(b1-3)GalNAc(a1-")
reduce_structure_level(glycan, to_level = "topological")
#> <glycan_structure[1]>
#> [1] Gal(??-?)GalNAc(??-
#> # Unique structures: 1
```
