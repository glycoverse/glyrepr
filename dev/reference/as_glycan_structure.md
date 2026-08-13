# Convert to Glycan Structure Vector

Convert an object to a glycan structure vector.

## Usage

``` r
as_glycan_structure(x, on_failure = c("error", "na"))
```

## Arguments

- x:

  An object to convert to a glycan structure vector. Can be an igraph
  object, a list of igraph objects, a character vector of
  IUPAC-condensed strings, or an existing glyrepr_structure object.

- on_failure:

  The failure policy for element-local parsing, validation, and
  canonicalization errors. `"error"` preserves the default strict
  behavior. `"na"` replaces failed elements with `NA` and emits one
  warning that reports their positions and failure reasons. Existing
  missing elements remain missing without a warning. Vector-level
  incompatibilities still produce an error.

## Value

A glyrepr_structure object.

## Details

Character input assumes the natural absolute configuration for
unprefixed monosaccharides. Less common configurations use a leading
`D-` or `L-`, such as `D-Fuc`, `L-Gul`, and `D-Fucf`.

Character input supports floating-part blocks before the main
IUPAC-condensed structure. `{Neu5Ac(a2-3)}<main>` allows every feasible
main-tree node as a candidate parent, while `{Neu5Ac(a2-3)|1,4}<main>`
restricts the candidates to main-tree nodes 1 and 4. Parent indices
follow residue order in the supplied main sequence and are remapped to
canonical IUPAC order in the result. The `|<parents>` suffix is a
`glyrepr` extension to curly-brace IUPAC notation. A singleton candidate
set is accepted as input but fully localizes the attachment, so
`{Neu5Ac(a2-3)|1}Gal(b1-4)GlcNAc(b1-` canonicalizes to the ordinary
structure `Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-`. An omitted parent list
behaves the same way when the main tree contains only one node.

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
as_glycan_structure(c("D-Fuc(a1-", "L-Gul(b1-", "D-Fucf(a1-"))
#> <glycan_structure[3]>
#> [1] D-Fuc(a1-
#> [2] L-Gul(b1-
#> [3] D-Fucf(a1-
#> # Unique structures: 3

# Parse a floating residue with two candidate parents
floating_iupac <- paste0(
  "{Neu5Ac(a2-3)|1,4}",
  "Gal(b1-4)GlcNAc(b1-2)Man(a1-3)",
  "[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]",
  "Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
)
as_glycan_structure(floating_iupac)
#> <glycan_structure[1]>
#> [1] {Neu5Ac(a2-3)|1,4}Gal(b1-4)GlcNAc(b1-2)Man(a1-3)[Gal(b1-4)GlcNAc(b1-2)Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> # Unique structures: 1

# Preserve valid elements while replacing an invalid element with NA
as_glycan_structure(
  c(valid = "Glc(?1-", invalid = "not-a-structure"),
  on_failure = "na"
)
#> Warning: 1 structure failed validation and was replaced with `NA`.
#> ✖ Position 2 (`invalid`): Could not parse IUPAC-condensed string:
#>   "not-a-structure" ℹ Invalid characters or format in IUPAC-condensed string
#> <glycan_structure[2]>
#> [1] valid    Glc(?1-
#> [2] invalid  NA
#> # Unique structures: 1
```
