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

## Data Structure Overview

A glycan structure vector is a vctrs vector with an additional S3 class
`glyrepr_structure`.

Each glycan structure must satisfy the following constraints:

### Graph Structure Requirements

- An ordinary structure must be a directed outward tree (reducing end as
  root).

- A structure with floating parts must be one annotated forest
  containing exactly one main outward tree and one outward tree per
  floating part.

- Floating substituents add graph metadata but no vertices or edges.

- Must have a graph attribute `anomer` in the format "a1" or "b1"

  - Unknown parts can be represented with "?", e.g., "?1", "a?", "??"

- May have a graph attribute `alditol`, containing one logical value.
  Missing attributes are treated as `FALSE` and canonicalized
  explicitly.

### Node Attributes

- `mono`: Monosaccharide names, must be known monosaccharide types

  - Generic names: Hex, HexNAc, dHex, NeuAc, etc.

  - Concrete names: Glc, Gal, Man, GlcNAc, etc.

  - Generic and concrete names may be mixed

  - NA values are not allowed

- `sub`: Substituent information

  - Single substituent format: "xY" (x = position, Y = substituent
    name), e.g., "2Ac", "3S"

  - Ambiguous substituent positions use slash-separated alternatives,
    e.g., "4/6S", "3/4/6Ac"

  - Multiple substituents separated by commas and ordered by position,
    e.g., "3Me,4Ac", "2S,6P"

  - Unknown substituent positions can be repeated, e.g., "?Me,?S"

  - No substituents represented by empty string ""

### Edge Attributes

- `linkage`: Glycosidic linkage information in format "a/bX-Y"

  - Standard format: e.g., "b1-4", "a2-3"

  - Unknown positions allowed: "a1-?", "b?-3", "??-?"

  - Partially unknown positions: "a1-3/6", "a1-3/6/9"

  - NA values are not allowed

### Floating Parts

Floating parts are disconnected substructures whose attachment to the
main tree is not fully localized. They are declared by the
`floating_parts` graph attribute, a list with one entry per floating
component. Each entry contains:

- `root`: the integer vertex index of the floating component root.

- `nodes`: all integer vertex indices in the floating component, ordered
  as the component appears in the complete IUPAC-condensed sequence.

- `linkage`: the virtual linkage from that root to its unresolved
  parent.

- `parents`: integer vertex indices outside the floating component. An
  empty integer vector means that all feasible nodes outside the
  component are candidates.

Canonical graphs always contain `nodes`. For backward compatibility,
input graphs may omit it; `glycan_structure()` derives the component
membership before validation and stores `nodes` in the canonical result.

During canonicalization, a floating part with exactly one effective
candidate parent is attached to that parent as an ordinary graph edge.
Attachments between floating components merge their `nodes` metadata and
can resolve further singleton domains. Only unresolved attachments
retain floating metadata, where the virtual attachment is metadata
rather than an edge and contributes to the canonical structure key.

### Floating Substituents

A floating substituent has known chemistry but an unresolved parent
residue. It is declared by the `floating_substituents` graph attribute,
a list with one entry per substituent. Each entry contains:

- `substituent`: one canonical substituent token such as `"6S"`,
  `"4/6Ac"`, or `"?Me"`.

- `parents`: integer residue vertex indices in the complete structure.
  An empty integer vector means that all feasible residue nodes are
  candidates.

A singleton candidate is normalized into the corresponding vertex's
`sub` attribute. Candidate parents must permit a conflict-free
assignment of occupied carbon positions. Floating-part assignments must
also be acyclic and connect every floating component to the main tree.

## Node and Edge Order

For an ordinary tree, the indices of vertices and linkages correspond
directly to their order in the printed IUPAC-condensed string. For
example, for the glycan
`Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-`, the vertices are
"Man", "Man", "Man", "GlcNAc", "GlcNAc", and the linkages are "a1-3",
"a1-6", "b1-4", "b1-4".

For a floating structure, floating-component vertices and edges precede
the main tree, exactly as their brace-enclosed components precede the
main glycan in the complete IUPAC-condensed string. Parent indices
written inside braces and stored in `floating_parts$parents` or
`floating_substituents$parents` use this same global order. Substituent
blocks contribute no vertex indices. A virtual floating attachment is
not an edge, and a floating substituent is not a vertex.

## NA Support

Glycan structure vectors support NA values for representing missing or
unknown structures:

- Create with `glycan_structure(NA)` or `glycan_structure(NULL)`

- Combine with valid structures: `c(struct1, NA, struct2)`

- Convert from character: `as_glycan_structure(c("Glc(a1-", NA))`

- `smap` functions skip NA elements gracefully

- [`is.na()`](https://rdrr.io/r/base/NA.html) returns `TRUE` for NA
  elements

## Naming Support

Glycan structure vectors can have names, which are preserved during
operations. This is particularly useful when working with the `glymotif`
package.

## Character conversion

A glycan structure vector is not a character vector. Use
[`as.character()`](https://rdrr.io/r/base/character.html) to explicitly
convert it to IUPAC-condensed strings when needed.

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

# Example 4: Parse a floating part with explicit candidate parents
floating <- as_glycan_structure(
  "{Neu5Ac(a2-3)|2,3}Gal(b1-3)[Gal(b1-4)]GlcNAc(a1-"
)
structure_floating_parts(floating)
#> # A tibble: 1 × 6
#>   glycan_id part_id root_node nodes     linkage parents  
#>       <int>   <int>     <int> <list>    <chr>   <list>   
#> 1         1       1         1 <int [1]> a2-3    <int [2]>

# Example 5: Parse a substituent with two candidate residues
floating_sub <- as_glycan_structure(
  "{6S|1,2}Gal(a1-3)Glc(a1-3)Man(a1-"
)
get_structure_graphs(floating_sub)$floating_substituents
#> [[1]]
#> [[1]]$substituent
#> [1] "6S"
#> 
#> [[1]]$parents
#> [1] 1 2
#> 
#> 

# Example 6: Check if object is a glycan structure
is_glycan_structure(simple_struct)  # TRUE
#> [1] TRUE
is_glycan_structure(graph)          # FALSE
#> [1] FALSE
```
