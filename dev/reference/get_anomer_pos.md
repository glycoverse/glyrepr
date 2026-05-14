# Get Anomer Positions

This function returns the anomer position for concrete monosaccharide
names.

## Usage

``` r
get_anomer_pos(mono)
```

## Arguments

- mono:

  A character vector of concrete monosaccharide names.

## Value

An integer vector of anomer positions.

## Examples

``` r
get_anomer_pos(c("Gal", "Neu5Ac"))
#> [1] 1 2
```
