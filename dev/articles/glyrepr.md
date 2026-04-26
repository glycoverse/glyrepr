# Getting Started with glyrepr

R gives us many data structures for different jobs. For example,
`c(1L, 2L, 3L)` represents an integer vector, `c("a", "b", "c")`
represents a character vector, and data frames represent tabular data.
Before we can inspect or manipulate a kind of data well, we first need a
suitable structure for representing it.

What about glycans?

When we talk about glycans, we usually talk about two related things:
their compositions and their structures. A composition can be
represented as a named vector, for example
`c(Hex = 3, HexNAc = 2, Neu5Ac = 1)`. A structure can be represented as
a graph, for example with the `igraph` package. Those representations
are useful, but managing the details by hand quickly becomes cumbersome,
especially for complex glycans or large datasets. This is where the
`glyrepr` package comes in.

`glyrepr` provides two main vector types:
[`glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_composition.md)
and
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md).
These vectors are designed to feel like ordinary R vectors: you can
subset them, concatenate them, sort them, put them in tibbles, and use
them in vectorized workflows.

These two representations are the foundation of `glycoverse`. Most
higher-level packages in the ecosystem build on them, so it is worth
taking a little time to get comfortable here.

``` r
library(glyrepr)
```

## Glycan Composition Vectors

Let’s start with glycan composition vectors.

### Creating Glycan Composition Vectors

A glycan composition vector can be created with either
[`glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_composition.md)
or
[`as_glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_composition.md).

[`glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_composition.md)
is the direct constructor. It takes one or more named vectors, where the
names are monosaccharides and the values are counts.

``` r
comps <- glycan_composition(
  c(Man = 5, GlcNAc = 2),
  c(Man = 3, Gal = 2, GlcNAc = 4),
  c(Man = 3, Gal = 2, GlcNAc = 4, Neu5Ac = 1, Fuc = 1)
)
comps
#> <glycan_composition[3]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

This creates a glycan composition vector called `comps` with three
glycan compositions.

[`as_glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_composition.md)
is more flexible and can convert several input types.

From a list of named vectors:

``` r
as_glycan_composition(list(
  c(Man = 5, GlcNAc = 2),
  c(Man = 3, Gal = 2, GlcNAc = 4),
  c(Man = 3, Gal = 2, GlcNAc = 4, Neu5Ac = 1, Fuc = 1)
))
#> <glycan_composition[3]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

This looks similar to
[`glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_composition.md),
but it is more convenient when your compositions are already stored in a
list, for example after reading or generating them programmatically.

``` r
comp_list <- list(
  c(Man = 5, GlcNAc = 2),
  c(Man = 3, Gal = 2, GlcNAc = 4),
  c(Man = 3, Gal = 2, GlcNAc = 4, Neu5Ac = 1, Fuc = 1)
)
as_glycan_composition(comp_list)
#> <glycan_composition[3]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

From a character vector:

``` r
as_glycan_composition(c("H5N2", "Hex(3)HexNAc(2)"))
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(3)HexNAc(2)
```

Both compact notation (`"H5N2"`) and Byonic-style notation
(`"Hex(3)HexNAc(2)"`) are supported.

From a glycan structure vector, which we will discuss below:

``` r
strucs <- c(o_glycan_core_1(), o_glycan_core_2())
as_glycan_composition(strucs)
#> <glycan_composition[2]>
#> [1] Gal(1)GalNAc(1)
#> [2] Gal(1)GlcNAc(1)GalNAc(1)
```

### Two Types of Monosaccharide Names

Before we move on, let’s briefly discuss the two types of monosaccharide
names used in `glyrepr`:

1.  **Generic names**: Hex, HexNAc, dHex, etc.
2.  **Specific names**: Man, GlcNAc, Fuc, etc.

Generic names are common in mass spectrometry data, where it is often
difficult to distinguish isomers from MS evidence alone. For example,
one `Hex` residue could be Man, Gal, Glc, or another hexose. Specific
names carry more biological detail and are common in glycan databases
and literature.

`glyrepr` supports both types of names, but you cannot mix them in the
same vector.

``` r
# This raises an error because the monosaccharide names are mixed.
try(as_glycan_composition(c("Hex(5)HexNAc(2)", "Man(5)GlcNAc(2)")), silent = TRUE)
```

The same rule applies to glycan structure vectors.

### Inspecting Glycan Composition Vectors

The main function for inspecting glycan compositions is
[`count_mono()`](https://glycoverse.github.io/glyrepr/dev/reference/count_mono.md).
Let’s demonstrate it using the `comps` vector we created earlier.

``` r
comps
#> <glycan_composition[3]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

The first argument is a glycan composition vector, and the second
argument is the monosaccharide you want to count.

``` r
count_mono(comps, "Man")
#> [1] 5 3 3
```

``` r
count_mono(comps, "Neu5Ac")
#> [1] 0 0 1
```

[`count_mono()`](https://glycoverse.github.io/glyrepr/dev/reference/count_mono.md)
works with both generic and specific monosaccharide names. It
understands the relationship between them, so it can count at the level
you ask for.

``` r
count_mono(comps, "Hex")
#> [1] 5 5 5
```

Note that both “Man” and “Gal” are counted as “Hex” here.

You can also omit the second argument to get the total monosaccharide
count for each composition.

``` r
count_mono(comps)
#> [1]  7  9 11
```

### Manipulating Glycan Composition Vectors

One useful mental model for glycan composition vectors (and glycan
structure vectors) is that they behave like atomic vectors. This means
they support familiar operations like subsetting, concatenation, and
sorting.

**Concatenation**:

``` r
c(comps, comps)
#> <glycan_composition[6]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
#> [4] Man(5)GlcNAc(2)
#> [5] Man(3)Gal(2)GlcNAc(4)
#> [6] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

**Subsetting**:

``` r
comps[1:2]
#> <glycan_composition[2]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
```

``` r
comps[integer()]
#> <glycan_composition[0]>
```

**Length**:

``` r
length(comps)
#> [1] 3
```

**Unique**:

``` r
dup_comps <- c(comps, comps)
dup_comps
#> <glycan_composition[6]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
#> [4] Man(5)GlcNAc(2)
#> [5] Man(3)Gal(2)GlcNAc(4)
#> [6] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

``` r
unique(dup_comps)
#> <glycan_composition[3]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

**Repeated**:

``` r
rep(comps, times = 2)
#> <glycan_composition[6]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
#> [4] Man(5)GlcNAc(2)
#> [5] Man(3)Gal(2)GlcNAc(4)
#> [6] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

**Sorting**:

``` r
sort(comps)
#> <glycan_composition[3]>
#> [1] Man(5)GlcNAc(2)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

``` r
sort(comps, decreasing = TRUE)
#> <glycan_composition[3]>
#> [1] Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
#> [2] Man(3)Gal(2)GlcNAc(4)
#> [3] Man(5)GlcNAc(2)
```

### Working with Tibbles

One of the most useful features of glycan composition vectors is that
they work smoothly with tibbles and data frames.

``` r
library(tibble)

tb <- tibble(
  id = c("glycan1", "glycan2", "glycan3"),
  composition = comps
)
tb
#> # A tibble: 3 × 2
#>   id      composition                         
#>   <chr>   <comp>                              
#> 1 glycan1 Man(5)GlcNAc(2)                     
#> 2 glycan2 Man(3)Gal(2)GlcNAc(4)               
#> 3 glycan3 Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)
```

You can use `tidyverse` functions to perform operations on the glycan
composition column.

``` r
library(dplyr)
#> 
#> Attaching package: 'dplyr'
#> The following objects are masked from 'package:stats':
#> 
#>     filter, lag
#> The following objects are masked from 'package:base':
#> 
#>     intersect, setdiff, setequal, union

tb |>
  mutate(n_sia = count_mono(composition, "Neu5Ac")) |>
  filter(n_sia > 0)
#> # A tibble: 1 × 3
#>   id      composition                          n_sia
#>   <chr>   <comp>                               <int>
#> 1 glycan3 Man(3)Gal(2)GlcNAc(4)Fuc(1)Neu5Ac(1)     1
```

### Missing Values and Names

Glycan composition vectors can also handle missing values.

``` r
comps_with_na <- glycan_composition(c(Man = 5, GlcNAc = 2), NA)
comps_with_na
#> <glycan_composition[2]>
#> [1] Man(5)GlcNAc(2)
#> [2] <NA>
```

``` r
count_mono(comps_with_na, "Man")
#> [1]  5 NA
```

For technical reasons, glycan composition vectors cannot have element
names right now. In practice, this is usually fine: the composition
itself remains the data value, and identifiers can live in a separate
tibble column.

## Glycan Structure Vectors

Now let’s move on to the core of glycan representation: glycan structure
vectors. Many of the same ideas apply, including the atomic vector
nature, vectorized operations, and seamless integration with tibbles.
Because structures are more complex than compositions, there are a few
additional features and considerations to keep in mind.

### Creating Glycan Structure Vectors

As with glycan composition vectors, you can create glycan structure
vectors with either
[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
or
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md).

Under the hood, glycan structures are represented as `igraph` objects,
because glycans are naturally graph-like. Most of the time, you do not
need to work with those graphs directly: `glyrepr` gives you a
higher-level interface for common structure operations.

[`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
is the direct constructor and takes one or more `igraph` objects.
Creating those graphs manually can be tedious, so you will usually start
with
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md)
instead. It can parse IUPAC-condensed strings into glycan structure
vectors.

If you are not familiar with IUPAC-condensed notation, check out this
[article](https://glycoverse.github.io/glyrepr/articles/iupac.html) for
a quick introduction. We recommend getting familiar with this notation,
because it is the main text language for communicating glycan structures
in `glycoverse`.

Parsing IUPAC-condensed strings is straightforward with
[`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md).
(For parsing other formats like GlycoCT, use the
[glyparse](https://github.com/glycoverse/glyparse) package.)

``` r
strucs <- as_glycan_structure(c(
  "Gal(b1-3)GalNAc(a1-",
  "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
))
strucs
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> # Unique structures: 2
```

### Inspecting Glycan Structure Vectors

`strucs` prints a bit like a character vector with colors, but under the
hood the glycan structures are stored as `igraph` objects. You can
access the underlying `igraph` objects using
[`get_structure_graphs()`](https://glycoverse.github.io/glyrepr/dev/reference/get_structure_graphs.md).

``` r
get_structure_graphs(strucs)
#> [[1]]
#> IGRAPH 81ecad3 DN-- 2 1 -- 
#> + attr: anomer (g/c), name (v/c), mono (v/c), sub (v/c), linkage (e/c)
#> + edge from 81ecad3 (vertex names):
#> [1] 2->1
#> 
#> [[2]]
#> IGRAPH 56a97d3 DN-- 3 2 -- 
#> + attr: anomer (g/c), name (v/c), mono (v/c), sub (v/c), linkage (e/c)
#> + edges from 56a97d3 (vertex names):
#> [1] 3->1 3->2
```

Details of how a glycan structure is modeled as a graph are covered in
the [glycan graph
vignette](https://glycoverse.github.io/glyrepr/articles/glycan-graph.html).
You do not need all of those details for everyday use, but a quick skim
can make the structure functions easier to understand.

The
[`count_mono()`](https://glycoverse.github.io/glyrepr/dev/reference/count_mono.md)
function also works with glycan structure vectors.

``` r
count_mono(strucs, "Gal")
#> [1] 1 1
```

Other functions, including
[`has_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/has_linkages.md),
[`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md),
and
[`get_structure_level()`](https://glycoverse.github.io/glyrepr/dev/reference/get_structure_level.md),
inspect specific aspects of the structures.

``` r
# This function works element-wise
has_linkages(strucs)
#> [1] TRUE TRUE
```

``` r
get_mono_type(strucs)
#> [1] "concrete"
```

``` r
get_structure_level(strucs)
#> [1] "intact"
```

### Structure Levels

The structures we have seen so far are “intact” structures. They contain
specific monosaccharides and complete linkage information. In real
datasets, though, glycan structures often have missing information. For
example, we might know the topology but not the linkages, or we might
only know generic monosaccharide classes.

To accommodate these scenarios, `glyrepr` defines four levels of glycan
structures:

- **Intact**: specific monosaccharides and complete linkage information.
- **Partial**: specific monosaccharides, with at least one missing
  linkage annotation.
- **Topological**: specific monosaccharides, but all linkage information
  is missing.
- **Basic**: generic monosaccharides, with linkage information treated
  as missing.

Structure levels are defined at the vector level, so one glycan
structure vector has one level.

Let’s see some examples.

**Intact structures**:

``` r
as_glycan_structure(c(
  "Gal(b1-3)GalNAc(a1-",
  "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
)) |> get_structure_level()
#> [1] "intact"
```

**Partial structures**:

``` r
as_glycan_structure(c(
  "Gal(b1-?)GalNAc(a1-",
  "Gal(b1-?)[GlcNAc(b1-6)]GalNAc(a1-"
)) |> get_structure_level()
#> [1] "partial"
```

**Topological structures**:

``` r
as_glycan_structure(c(
  "Gal(??-?)GalNAc(??-",
  "Gal(??-?)[GlcNAc(??-?)]GalNAc(??-"
)) |> get_structure_level()
#> [1] "topological"
```

**Basic structures**:

``` r
as_glycan_structure(c(
  "Hex(??-?)HexNAc(??-",
  "Hex(??-?)[HexNAc(??-?)]HexNAc(??-"
)) |> get_structure_level()
#> [1] "basic"
```

In theory, you can have something like `"Hex(b1-3)HexNAc(a1-"`, with
generic monosaccharides but all linkages intact. In practice, linkage
information is usually harder to obtain than monosaccharide identities,
so this situation is rare.

If you create such a vector, `glyrepr` classifies it as `"basic"` and
warns you.

``` r
as_glycan_structure(c(
  "Hex(a1-3)HexNAc(a1-",
  "Hex(a1-3)[HexNAc(b1-6)]HexNAc(a1-"
)) |> get_structure_level()
#> Warning: Generic glycan structures with linkage annotations are treated as "basic".
#> ℹ Linkage information is ignored when residues are generic.
#> [1] "basic"
```

### Manipulating Glycan Structure Vectors

Glycan structure vectors also support vectorized operations like
subsetting and concatenation. They also have structure-specific helpers
for reducing resolution, removing linkages, and removing substituents.

``` r
strucs
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> # Unique structures: 2
```

``` r
reduce_structure_level(strucs, to_level = "basic")
#> <glycan_structure[2]>
#> [1] Hex(??-?)HexNAc(??-
#> [2] HexNAc(??-?)[Hex(??-?)]HexNAc(??-
#> # Unique structures: 2
```

``` r
remove_linkages(strucs)  # same as reduce_structure_level(strucs, to_level = "topological")
#> <glycan_structure[2]>
#> [1] Gal(??-?)GalNAc(??-
#> [2] GlcNAc(??-?)[Gal(??-?)]GalNAc(??-
#> # Unique structures: 2
```

``` r
strucs_with_subs <- as_glycan_structure(c(
  "Gal6S(b1-3)GalNAc(a1-",
  "Gal6S(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
))
remove_substituents(strucs_with_subs)
#> <glycan_structure[2]>
#> [1] Gal(b1-3)GalNAc(a1-
#> [2] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> # Unique structures: 2
```

### Missing Values and Names

Like glycan composition vectors, glycan structure vectors can also
handle missing values.

Different from glycan composition vectors, though, glycan structure
vectors can have names. This is a very useful feature we introduced in
version 0.10.0.

``` r
strings <- c(
  glycan1 = "Gal(b1-3)GalNAc(a1-",
  glycan2 = "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"
)
as_glycan_structure(strings)
#> <glycan_structure[2]>
#> [1] glycan1  Gal(b1-3)GalNAc(a1-
#> [2] glycan2  Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> # Unique structures: 2
```

You will find the names useful when working with the `glymotif` package.

## What’s Next?

This vignette sets the foundation for working with glycan
representations in `glycoverse`. From here, you can explore the packages
that build on these representations:

- [glyparse](https://github.com/glycoverse/glyparse): parsing text
  nomenclatures into glycan representations.
- [glymotif](https://github.com/glycoverse/glymotif): identifying and
  analyzing glycan motifs.
- [glydet](https://github.com/glycoverse/glydet): calculating derived
  traits and quantifying motifs.
- [glydraw](https://github.com/glycoverse/glydraw): visualizing glycan
  structures with SNFG notation.
- [glyenzy](https://github.com/glycoverse/glyenzy): inspecting and
  modeling glycan biosynthesis.

Happy exploring.
