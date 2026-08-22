# Get Alditol Status

Determine whether the reducing-end residue of a glycan is an alditol.
Alditol status is stored as the graph-level logical attribute `alditol`.
Legacy graphs without this attribute are treated as non-alditols.

## Usage

``` r
get_alditol(x)
```

## Arguments

- x:

  A glycan structure vector or one glycan `igraph`.

## Value

A logical vector for structure-vector input, or one logical scalar for
graph input. Missing structure-vector elements return `NA`.

## Examples

``` r
glycans <- as_glycan_structure(c(
  reduced = "Gal(b1-4)GlcNAc-ol(a1-",
  ordinary = "Gal(b1-4)GlcNAc(a1-"
))
get_alditol(glycans)
#>  reduced ordinary 
#>     TRUE    FALSE 
```
