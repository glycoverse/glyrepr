# Remove All Substituents from a Glycan

This function replaces all substituents in a glycan structure with empty
strings.

## Usage

``` r
remove_substituents(glycan)
```

## Arguments

- glycan:

  A glyrepr_structure vector.

## Value

A glyrepr_structure vector with all substituents removed.

## Examples

``` r
(glycan <- glycan_structure(o_glycan_core_1()))
#> <glycan_structure[1]>
#> [1] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
remove_substituents(glycan)
#> <glycan_structure[1]>
#> [1] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 1
```
