# Infer Anomer Positions

This function infers the anomer position for concrete or generic
monosaccharide names.

## Usage

``` r
infer_anomer_pos(mono)

get_anomer_pos(mono)
```

## Arguments

- mono:

  A character vector of monosaccharide names.

## Value

An integer vector of anomer positions.

## Examples

``` r
infer_anomer_pos(c("Gal", "Hex", "Neu5Ac"))
#> [1] 1 1 2
```
