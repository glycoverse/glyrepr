# Generate Possible Linkages

Given an obscure linkage format (having "?", e.g. "a2-?"), this function
generates all possible linkages based on the format. See
[`valid_linkages()`](https://glycoverse.github.io/glyrepr/reference/valid_linkages.md)
for details.

The ranges of possible anomers, first positions, and second positions
can be specified using `anomer_range`, `pos1_range`, and `pos2_range`.

## Usage

``` r
possible_linkages(
  linkage,
  anomer_range = c("a", "b"),
  pos1_range = 1:2,
  pos2_range = 1:9,
  include_unknown = FALSE
)
```

## Arguments

- linkage:

  A linkage string.

- anomer_range:

  A character vector of possible anomers. Default is `c("a", "b")`.

- pos1_range:

  A numeric vector of possible first positions. Default is `1:2`.

- pos2_range:

  A numeric vector of possible second positions. Default is `1:9`.

- include_unknown:

  A logical value. If `TRUE`, "?" will be included. Default is `FALSE`.

## Value

A character vector of possible linkages.

## See also

[`has_linkages()`](https://glycoverse.github.io/glyrepr/reference/has_linkages.md),
[`remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.md),
[`valid_linkages()`](https://glycoverse.github.io/glyrepr/reference/valid_linkages.md)

## Examples

``` r
possible_linkages("a2-?")
#> [1] "a2-1" "a2-2" "a2-3" "a2-4" "a2-5" "a2-6" "a2-7" "a2-8" "a2-9"
possible_linkages("??-2")
#> [1] "a1-2" "b1-2" "a2-2" "b2-2"
possible_linkages("a1-3")
#> [1] "a1-3"
possible_linkages("a?-?", pos1_range = 2, pos2_range = c(2, 3))
#> [1] "a2-2" "a2-3"
possible_linkages("?1-6", include_unknown = TRUE)
#> [1] "a1-6" "b1-6" "?1-6"
```
