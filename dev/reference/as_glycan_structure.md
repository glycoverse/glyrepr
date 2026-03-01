# Convert to Glycan Structure Vector

Convert an object to a glycan structure vector.

## Usage

``` r
as_glycan_structure(x)
```

## Arguments

- x:

  An object to convert to a glycan structure vector. Can be an igraph
  object, a list of igraph objects, a character vector of
  IUPAC-condensed strings, or an existing glyrepr_structure object.

## Value

A glyrepr_structure object.

## Examples

``` r
library(igraph)
#> 
#> Attaching package: ‘igraph’
#> The following objects are masked from ‘package:stats’:
#> 
#>     decompose, spectrum
#> The following object is masked from ‘package:base’:
#> 
#>     union

# Convert a single igraph
graph <- make_graph(~ 1-+2)
V(graph)$mono <- c("GlcNAc", "GlcNAc")
V(graph)$sub <- ""
E(graph)$linkage <- "b1-4"
graph$anomer <- "a1"
as_glycan_structure(graph)
#> <glycan_structure[1]>
#> [1] GlcNAc(b1-4)GlcNAc(a1-
#> # Unique structures: 1

# Convert a list of igraphs
o_glycan_vec <- o_glycan_core_1()
o_glycan_graph <- get_structure_graphs(o_glycan_vec)
as_glycan_structure(list(graph, o_glycan_graph))
#> <glycan_structure[2]>
#> [1] GlcNAc(b1-4)GlcNAc(a1-
#> [2] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 2

# Convert a character vector of IUPAC-condensed strings
as_glycan_structure(c("GlcNAc(b1-4)GlcNAc(b1-", "Man(a1-2)GlcNAc(b1-"))
#> <glycan_structure[2]>
#> [1] GlcNAc(b1-4)GlcNAc(b1-
#> [2] Man(a1-2)GlcNAc(b1-
#> # Unique structures: 2
```
