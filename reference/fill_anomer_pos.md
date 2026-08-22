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
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  vector or glycan `igraph` with concrete or generic monosaccharides.

## Value

An object of the same representation as `strucs` with anomer positions
added where missing. Graph input retains its vertex IDs and order.

## Details

For anomer positions that are already specified in the input structures,
this function does not modify them.

For a structure with floating parts, the reducing-end position is
inferred from the root of the main tree. Positions in virtual
floating-part attachment linkages are inferred from each floating part's
root residue.

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
