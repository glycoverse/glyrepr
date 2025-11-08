# Determine if a Glycan Structure has Linkages

Unknown linkages in a glycan structure are represented by "??-?". Also,
a linkage can be partially known (e.g. "a?-?"). This function checks if
a glycan structure has linkages, in a strict or lenient way.

## Usage

``` r
has_linkages(glycan, strict = FALSE)
```

## Arguments

- glycan:

  A
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  vector.

- strict:

  A logical value.

  - If `FALSE` (default), a glycan is considered to have linkages if any
    linkage is partially known (not "??-?").

  - If `TRUE`, a glycan is considered to have linkages only if all
    linkages are fully determined (no "?" in the linkage).

## Value

A logical vector indicating if each glycan structure has linkages.

## See also

[`remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.md),
[`possible_linkages()`](https://glycoverse.github.io/glyrepr/reference/possible_linkages.md)

## Examples

``` r
glycan <- o_glycan_core_1(linkage = TRUE)
has_linkages(glycan)
#> [1] TRUE
print(glycan)
#> <glycan_structure[1]>
#> [1] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1

glycan <- remove_linkages(glycan)
has_linkages(glycan)
#> [1] FALSE
print(glycan)
#> <glycan_structure[1]>
#> [1] Gal(??-?)GalNAc(??-
#> # Unique structures: 1

glycan <- as_glycan_structure("Gal(b1-?)GalNAc(a1-")
has_linkages(glycan)
#> [1] TRUE
has_linkages(glycan, strict = TRUE)
#> [1] FALSE
```
