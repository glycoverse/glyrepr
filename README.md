
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyrepr

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/glyrepr)](https://CRAN.R-project.org/package=glyrepr)
[![R-CMD-check](https://github.com/fubin1999/glyrepr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/fubin1999/glyrepr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/fubin1999/glyrepr/graph/badge.svg)](https://app.codecov.io/gh/fubin1999/glyrepr)
<!-- badges: end -->

This package is a wrapper of ‘igraph’ for representing glycan structures
in R. It provides a set of functions to create and manipulate glycan
structures, and to visualize them in a graph format.

## Installation

You can install the development version of glyrepr from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("fubin1999/glyrepr")
```

## Example

``` r
library(glyrepr)

glycan <- n_glycan_core()
print(glycan, verbose = TRUE)
#> Glycan Graph (NE)
#> GlcNAc: 2, Man: 3
#> ------------------
#> GlcNAc
#> └─GlcNAc (b1-4)
#>   └─Man (b1-4)
#>     ├─Man (a1-3)
#>     └─Man (a1-6)
```
