# Create a Glycan Structure Vector

`glycan_structure()` creates an efficient glycan structure vector for
storing and processing glycan molecular structures. The function employs
hash-based deduplication mechanisms, making it suitable for
glycoproteomics, glycomics analysis, and glycan structure comparison
studies.

## Usage

``` r
glycan_structure(...)

is_glycan_structure(x)
```

## Arguments

- ...:

  igraph graph objects to be converted to glycan structures, or existing
  glycan structure vectors. Supports mixed input of multiple objects.

- x:

  An object to check or convert.

## Value

A `glyrepr_structure` class glycan structure vector object.

## Core Features

- **Efficient Storage**: Uses hash values of IUPAC codes for
  deduplication, avoiding redundant storage of identical glycan
  structures

- **Graph Model Representation**: Each glycan structure is represented
  as a directed graph where nodes are monosaccharides and edges are
  glycosidic linkages

- **Vectorized Operations**: Supports R's vectorized operations for
  batch processing of glycan data

- **Type Safety**: Built on the vctrs package, providing type-safe
  operations

## Data Structure Overview

A glycan structure vector is a vctrs record with an additional S3 class
`glyrepr_structure`. Therefore, `sloop::s3_class()` returns the class
hierarchy `c("glyrepr_structure", "vctrs_rcrd")`.

Each glycan structure must satisfy the following constraints:

### Graph Structure Requirements

- Must be a directed graph with an outward tree structure (reducing end
  as root)

- Must have a graph attribute `anomer` in the format "a1" or "b1"

  - Unknown parts can be represented with "?", e.g., "?1", "a?", "??"

### Node Attributes

- `mono`: Monosaccharide names, must be known monosaccharide types

  - Generic names: Hex, HexNAc, dHex, NeuAc, etc.

  - Concrete names: Glc, Gal, Man, GlcNAc, etc.

  - Cannot mix generic and concrete names

  - NA values are not allowed

- `sub`: Substituent information

  - Single substituent format: "xY" (x = position, Y = substituent
    name), e.g., "2Ac", "3S"

  - Multiple substituents separated by commas and ordered by position,
    e.g., "3Me,4Ac", "2S,6P"

  - No substituents represented by empty string ""

### Edge Attributes

- `linkage`: Glycosidic linkage information in format "a/bX-Y"

  - Standard format: e.g., "b1-4", "a2-3"

  - Unknown positions allowed: "a1-?", "b?-3", "??-?"

  - Partially unknown positions: "a1-3/6", "a1-3/6/9"

  - NA values are not allowed

## Node and Edge Order

The indices of vertices and linkages in a glycan correspond directly to
their order in the IUPAC-condensed string, which is printed when you
print a `glycan_structure()`. For example, for the glycan
`Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-`, the vertices are
"Man", "Man", "Man", "GlcNAc", "GlcNAc", and the linkages are "a1-3",
"a1-6", "b1-4", "b1-4".

## Use Cases

- **Glycoproteomics Analysis**: Processing glycan structure information
  from mass spectrometry data

- **Glycomics Research**: Comparing glycan expression profiles across
  different samples or conditions

- **Structure-Function Analysis**: Studying relationships between glycan
  structures and biological functions

- **Database Queries**: Performing structure matching and searches in
  glycan databases

## Examples

``` r
library(igraph)

# Example 1: Create a simple glycan structure GlcNAc(b1-4)GlcNAc
graph <- make_graph(~ 1-+2)  # Create graph with two monosaccharides
V(graph)$mono <- c("GlcNAc", "GlcNAc")  # Set monosaccharide types
V(graph)$sub <- ""  # No substituents
E(graph)$linkage <- "b1-4"  # b1-4 glycosidic linkage
graph$anomer <- "a1"  # a anomeric carbon

# Create glycan structure vector
simple_struct <- glycan_structure(graph)
print(simple_struct)
#> <glycan_structure[1]>
#> [1] GlcNAc(b1-4)GlcNAc(a1-
#> # Unique structures: 1

# Example 2: Use predefined glycan core structures
n_core <- n_glycan_core()  # N-glycan core structure
o_core1 <- o_glycan_core_1()  # O-glycan Core 1 structure

# Create vector with multiple structures
multi_struct <- glycan_structure(n_core, o_core1)
print(multi_struct)
#> <glycan_structure[2]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [2] Gal(b1-3)GalNAc(a1-
#> # Unique structures: 2

# Example 3: Create complex structure with substituents
complex_graph <- make_graph(~ 1-+2-+3)
V(complex_graph)$mono <- c("GlcNAc", "Gal", "Neu5Ac")
V(complex_graph)$sub <- c("", "", "")  # Add substituents as needed
E(complex_graph)$linkage <- c("b1-4", "a2-3")
complex_graph$anomer <- "b1"

complex_struct <- glycan_structure(complex_graph)
print(complex_struct)
#> <glycan_structure[1]>
#> [1] Neu5Ac(a2-3)Gal(b1-4)GlcNAc(b1-
#> # Unique structures: 1

# Example 4: Check if object is a glycan structure
is_glycan_structure(simple_struct)  # TRUE
#> [1] TRUE
is_glycan_structure(graph)          # FALSE
#> [1] FALSE

# Example 5: Mix different input types
mixed_struct <- glycan_structure(graph, o_glycan_core_2(), simple_struct)
print(mixed_struct)
#> <glycan_structure[3]>
#> [1] GlcNAc(b1-4)GlcNAc(a1-
#> [2] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> [3] GlcNAc(b1-4)GlcNAc(a1-
#> # Unique structures: 2
```
