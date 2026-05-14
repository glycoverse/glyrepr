# Fill Anomer Positions

Add anomer positions to glycan structures with missing anomer position
information. For example, `"Gal(??-?)GalNAc(??-"` is converted to
`"Gal(?1-?)GalNAc(?1-"`.

## Usage

``` r
fill_anomer_pos(strucs)
```

## Arguments

- strucs:

  A
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
  vector with concrete monosaccharides.

## Value

A
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
vector with anomer positions added where missing.

## Details

For anomer positions that are already specified in the input structures,
this function does not modify them.

## Examples

``` r
glycans <- as_glycan_structure(c(
  "Gal(??-?)GalNAc(??-",
  "Neu5Ac(??-?)Gal(??-?)GalNAc(??-"
))
fill_anomer_pos(glycans)
#> <glycan_structure[2]>
#> [1] Gal(?1-?)GalNAc(?1-
#> [2] Neu5Ac(?2-?)Gal(?1-?)GalNAc(?1-
#> # Unique structures: 2
```
