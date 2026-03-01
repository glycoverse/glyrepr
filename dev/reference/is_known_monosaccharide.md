# Check if a Monosaccharide is Known

This function checks if a vector of monosaccharide names are known.

## Usage

``` r
is_known_monosaccharide(mono)
```

## Arguments

- mono:

  A character vector of monosaccharide names.

## Value

A logical vector.

## Examples

``` r
is_known_monosaccharide(c("Gal", "Hex"))
#> [1] TRUE TRUE
is_known_monosaccharide(c("X", "Hx", "Nac"))
#> [1] FALSE FALSE FALSE
```
