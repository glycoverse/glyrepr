# Get the Anomeric information

Get the Anomeric information

## Usage

``` r
get_anomer(x)
```

## Arguments

- x:

  A glycan structure vector or a glycan `igraph`.

## Value

A character vector for structure-vector input, or a character scalar for
graph input.

## Examples

``` r
x <- n_glycan_core()
get_anomer(x)
#> [1] "b1"
```
