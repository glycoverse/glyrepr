# Remove All Linkages from a Glycan

This function replaces all linkages in a glycan structure with "??-?",
as well as the reducing end anomer with "??-".

## Usage

``` r
remove_linkages(glycan)
```

## Arguments

- glycan:

  A glyrepr_structure vector.

## Value

A glyrepr_structure vector with all linkages removed.

## Examples

``` r
glycan <- o_glycan_core_1(linkage = TRUE)
glycan
#> <glycan_structure[1]>
#> [1] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
remove_linkages(glycan)
#> <glycan_structure[1]>
#> [1] Gal(??-?)GalNAc(??-
#> # Unique structures: 1
```
