
<!-- README.md is generated from README.Rmd. Please edit that file -->

# glyrepr <a href="https://glycoverse.github.io/glyrepr/"><img src="man/figures/logo.png" align="right" height="138" /></a>

<!-- badges: start -->

[![Lifecycle:
experimental](https://img.shields.io/badge/lifecycle-experimental-orange.svg)](https://lifecycle.r-lib.org/articles/stages.html#experimental)
[![CRAN
status](https://www.r-pkg.org/badges/version/glyrepr)](https://CRAN.R-project.org/package=glyrepr)
[![R-CMD-check](https://github.com/glycoverse/glyrepr/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/glycoverse/glyrepr/actions/workflows/R-CMD-check.yaml)
[![Codecov test
coverage](https://codecov.io/gh/glycoverse/glyrepr/graph/badge.svg)](https://app.codecov.io/gh/glycoverse/glyrepr)
<!-- badges: end -->

This package is a wrapper of ‘igraph’ for representing glycan structures
in R. It introduces an S3 class called “glycan_graph” along with a suite
of functions to create and manipulate glycan structures.

If you are a glycomics or glycoproteomics researcher, you likely won’t
need to interact with this package directly. All you need to know is
that it defines the “glycan_graph” S3 class to represent glycan
structures, so you won’t be surprised if you encounter this name
elsewhere in `glycoverse.`

In fact, the functions in this packages are heavily used by other
`glycoverse` packages such as
[glyparse](https://github.com/glycoverse/glyparse) and
[glymotif](https://github.com/glycoverse/glymotif), which you are
probably already using or will use.

## Installation

You can install the development version of glyrepr from
[GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("glycoverse/glyrepr")
```

## Example

``` r
library(glyrepr)

glycan <- n_glycan_core()
print(glycan, verbose = TRUE)
#> Glycan Graph (NE)
#> GlcNAc: 2, Man: 3
#> ------------------
#> GlcNAc (?1-)
#> └─GlcNAc (b1-4)
#>   └─Man (b1-4)
#>     ├─Man (a1-3)
#>     └─Man (a1-6)
```
