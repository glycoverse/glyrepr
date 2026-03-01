# Check if Linkages are Valid

Valid linkages are in the form of "a1-2", "b1-4", "a?-1", etc.
Specifically, the pattern is `xy-z`:

- `x`: the anomer, either "a", "b", or "?".

- `y`: the first position, either "1", "2" or "?".

- `z`: the second position, either a 1-9 digit or "?". Can also be
  multiple positions separated by "/", e.g. "1/2/3". "?" could not be
  used with "/".

## Usage

``` r
valid_linkages(linkages)
```

## Arguments

- linkages:

  A character vector of linkages.

## Value

A logical vector.

## Examples

``` r
# Valid linkages
valid_linkages(c("a1-2", "?1-4", "a?-1", "b?-?", "??-?", "a1/2-3"))
#> [1]  TRUE  TRUE  TRUE  TRUE  TRUE FALSE

# Invalid linkages
valid_linkages(c("a1-2/?", "1-4", "a/b1-2", "c1-2", "a9-1"))
#> [1] FALSE FALSE FALSE FALSE FALSE
```
