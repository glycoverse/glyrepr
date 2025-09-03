
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

`glyrepr` provides two important representations of glycans: glycan
compositions and glycan structures. This package is the core of
`glycoverse` ecosystem, as it provides the basic data structures and
functions for representing and manipulating glycans.

In fact, the functions in this packages are heavily used by other
`glycoverse` packages such as
[glyparse](https://github.com/glycoverse/glyparse) and
[glymotif](https://github.com/glycoverse/glymotif), which you are
probably already using or will use.

## Installation

You can install the latest release of glyrepr from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("glycoverse/glyrepr@*release")
```

Or install the development version:

``` r
devtools::install_github("glycoverse/glyrepr")
```

## Documentation

-   🚀 Get started:
    [Here](https://glycoverse.github.io/glyrepr/articles/glyrepr.html)
-   🔍 Underlying graph representation:
    [Here](https://glycoverse.github.io/glyrepr/articles/glycan-graph.html)
-   🔧 Faster structure operations:
    [Here](https://glycoverse.github.io/glyrepr/articles/smap.html)
-   ✏️ IUPAC-condensed glycan text representation:
    [Here](https://glycoverse.github.io/glyrepr/articles/iupac.html)
-   📚 Reference:
    [Here](https://glycoverse.github.io/glyrepr/reference/index.html)

## Role in `glycoverse`

As the cornerstone of the `glycoverse` ecosystem, this package provides
two fundamental data structures for representing glycans:
`glycan_composition()` and `glycan_structure()`. These specialized data
types are what distinguish `glycoverse` from other omics analysis
frameworks, serving as the foundation for higher-level packages like
`glymotif`, which build upon them to perform advanced glycan analysis.

## Example

``` r
library(glyrepr)

# Create glycan compositions
glycan_composition(
  c(Man = 5, GlcNAc = 2),
  c(Man = 3, Gal = 2, GlcNAc = 4, Fuc = 1, Neu5Ac = 2)
)
#> <glycan_composition[2]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(2)

# Parse IUPAC-condensed glycan text representation
# `glyrepr` supports IUPAC-condensed glycan text representation natively.
# For other formats like WURCS or glycoCT, use the `glyparse` package.
# For example, the following two glycan structures are equivalent:
structures <- as_glycan_structure(c("GlcNAc(b1-4)GlcNAc(?1-", "Man(a1-2)GlcNAc(?1-"))

# Get the composition of a glycan structure
as_glycan_composition(structures)
#> <glycan_composition[2]>
#> [1] GlcNAc(2)
#> [2] Man(1)GlcNAc(1)

# Count the number of given residues
count_mono(structures, "Hex")
#> [1] 0 1
count_mono(glycan_composition(c(Man = 3, GlcNAc = 2, Gal = 2)), "Hex")
#> [1] 5
```
