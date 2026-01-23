# Getting Started with glyrepr

Welcome to the world of glycan analysis! If you’ve ever tried to work
with glycans computationally, you know the struggle: these tree-like
molecules are notoriously difficult to represent and analyze compared to
their linear cousins like proteins or DNA. That’s where `glyrepr` comes
to the rescue.

Think of `glyrepr` as your **glycan translator** — it teaches your
computer how to “speak glycan” fluently. Whether you’re dealing with
compositions (what’s in the glycan) or structures (how it’s connected),
this package has got you covered.

``` r
library(glyrepr)
```

## Quick Start: What Are We Talking About?

Before we dive in, let’s establish our vocabulary. Don’t worry — it’s
simpler than it sounds!

| Term               | What It Means                                   | Example                    |
|--------------------|-------------------------------------------------|----------------------------|
| **Composition**    | The “ingredients list” — how many of each sugar | `Hex(5)HexNAc(2)`          |
| **Structure**      | The “blueprint” — how sugars are connected      | `Man(a1-3)Man(b1-4)GlcNAc` |
| **Monosaccharide** | A single sugar unit (the building blocks)       | `Gal`, `Man`, `Hex`        |
| **Linkage**        | The “glue” between sugars                       | `a1-3`, `b1-4`             |
| **Substitution**   | Chemical decorations on sugars                  | `6Ac`, `3Me`               |

🔍 **Pro tip**: We distinguish between **generic** sugars (like mystery
boxes labeled “Hex”) and **concrete** sugars (like specific boxes
labeled “Galactose”).

## Part 1: Compositions — The Easy Start

Let’s start with something straightforward: glycan compositions. Think
of these as ingredient lists for your favorite recipes.

### Creating Your First Compositions

There are three ways to create compositions, each with its own
superpower:

**Method 1: The Direct Approach**

``` r
# Just tell R what you have
glycan_composition(c(Man = 5, GlcNAc = 2), c(Gal = 1, GalNAc = 1))
#> <glycan_composition[2]>
#> [1] Man(5)GlcNAc(2)
#> [2] Gal(1)GalNAc(1)
```

**Method 2: The Programmatic Way**

``` r
# Perfect when you're processing data from files or databases
comp_list <- list(c(Man = 5, GlcNAc = 2), c(Gal = 1, GalNAc = 1))
as_glycan_composition(comp_list)
#> <glycan_composition[2]>
#> [1] Man(5)GlcNAc(2)
#> [2] Gal(1)GalNAc(1)
```

**Method 3: The Parser**

``` r
# Copy-paste from your mass spec software? No problem!
as_glycan_composition(c("Hex(5)HexNAc(2)", "H1N1"))
#> <glycan_composition[2]>
#> [1] Hex(5)HexNAc(2)
#> [2] Hex(1)HexNAc(1)
```

### The Magic of Colors 🌈

Here’s something cool: when you run these examples in your R console,
you’ll see the concrete monosaccharides (like Gal and GalNAc) displayed
in beautiful colors! These follow the [SNFG
standard](https://www.ncbi.nlm.nih.gov/glycans/snfg.html) — the
universal “color code” for glycans. Think of it as the glycan rainbow
🌈.

### Smart Counting with `count_mono()`

Now here’s where `glyrepr` shows its intelligence:

``` r
comp <- glycan_composition(c(Gal = 1, Man = 1, GalNAc = 1))

# How many galactose residues?
count_mono(comp, "Gal")
#> [1] 1

# How many hexose residues? (This includes Gal and Man!)
count_mono(comp, "Hex")
#> [1] 2
```

Notice how
[`count_mono()`](https://glycoverse.github.io/glyrepr/reference/count_mono.md)
is smart enough to know that galactose and mannose are both hexoses?
That’s the power of understanding glycan hierarchies!

## Part 2: Structures — Where the Magic Happens

Compositions are nice, but structures are where `glyrepr` truly shines.
This is like going from knowing the ingredients to understanding the
actual recipe and cooking method.

### Your First Glycan Structures

Let’s work with some real glycan structures. These strings below are
called the “IUPAC-condensed” glycan text representations. They might
look cryptic, but they’re actually quite readable once you get the hang
of it. To learn about them, check out [this
article](https://glycoverse.github.io/glyrepr/articles/iupac.html).

``` r
iupacs <- c(
  "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-",  # The famous N-glycan core
  "Gal(b1-3)GalNAc(a1-",                                  # O-glycan core 1
  "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-",                    # O-glycan core 2
  "Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-",      # A branched mannose tree
  "GlcNAc6Ac(b1-4)Glc3Me(a1-"                             # With some decorations
)

struc <- as_glycan_structure(iupacs)
struc
#> <glycan_structure[5]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [2] Gal(b1-3)GalNAc(a1-
#> [3] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> [4] Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-
#> [5] GlcNAc6Ac(b1-4)Glc3Me(a1-
#> # Unique structures: 5
```

### The Secret Sauce: Unique Structure Optimization

Here’s where `glyrepr` gets really clever. Notice that “# Unique
structures: 5” message? This isn’t just informational — it’s the key to
lightning-fast performance.

Let’s see this optimization in action:

``` r
# Create a big dataset with lots of repetition
large_struc <- rep(struc, 1000)  # 5,000 structures total
large_struc
#> <glycan_structure[5000]>
#> [1] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [2] Gal(b1-3)GalNAc(a1-
#> [3] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> [4] Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-
#> [5] GlcNAc6Ac(b1-4)Glc3Me(a1-
#> [6] Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-
#> [7] Gal(b1-3)GalNAc(a1-
#> [8] Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-
#> [9] Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-
#> [10] GlcNAc6Ac(b1-4)Glc3Me(a1-
#> ... (4990 more not shown)
#> # Unique structures: 5
```

Still showing “# Unique structures: 5”! This means `glyrepr` is storing
only 5 unique graphs internally, not 5,000. This is like having a smart
library system that stores only one copy of each book, no matter how
many people want to read it.

### Performance That Will Blow Your Mind 🚀

Let’s put this to the test:

``` r
library(tictoc)

tic("Converting 5 structures")
result_small <- convert_to_generic(struc)
toc()
#> Converting 5 structures: 0.022 sec elapsed

tic("Converting 5,000 structures")
result_large <- convert_to_generic(large_struc)
toc()
#> Converting 5,000 structures: 0.021 sec elapsed
```

**Mind = blown!** 🤯 The performance is nearly identical because
`glyrepr` only processes each unique structure once, then cleverly
expands the results.

### Understanding Structure Resolution Levels 🔬

Not all glycan structures are created equal — they come in different
levels of detail, like zoom levels on a map. `glyrepr` recognizes four
resolution levels:

- **“intact”**: The full picture — all monosaccharides are concrete
  (e.g., “Man”, “GlcNAc”), and all linkages are fully determined (e.g.,
  “a2-3”, “b1-4”).
- **“partial”**: Almost there — all monosaccharides are concrete (e.g.,
  “Man”, “GlcNAc”), but some linkage information is missing (e.g.,
  “a2-?”).
- **“topological”**: We know what’s there, but not how they connect —
  all monosaccharides are concrete (e.g., “Man”, “GlcNAc”), but the
  linkage information is completely unknown (“??-?”).
- **“basic”**: The minimalist view — all monosaccharides are generic
  (e.g., “Hex”, “HexNAc”), and the linkage information is completely
  unknown (“??-?”).

💡 **Fun fact**: In theory, you could have a glycan with generic
monosaccharides but fully determined linkages (e.g.,
“Hex(b1-3)HexNAc(a1-”)). But in practice, this is almost unheard of —
linkage information is much harder to obtain than monosaccharide
information. That’s why `glyrepr` assigns these structures to the
“basic” level too.

You can get the structure level for a glycan structure vector with
[`get_structure_level()`](https://glycoverse.github.io/glyrepr/reference/get_structure_level.md):

``` r
# Concrete structures (various linkage detail levels)
concrete_glycans <- as_glycan_structure(c(
  "Gal(b1-3)GalNAc(a1-",
  "Gal(b1-?)GalNAc(a1-",
  "Gal(??-?)GalNAc(??-"
))
get_structure_level(concrete_glycans)
#> [1] "intact"      "partial"     "topological"

# Generic structures
generic_glycans <- as_glycan_structure(c(
  "Hex(??-?)HexNAc(??-",
  "Hex(b1-3)HexNAc(a1-"
))
get_structure_level(generic_glycans)
#> [1] "basic" "basic"
```

### Structure Manipulation Tools

`glyrepr` comes with several handy tools for structure manipulation:

**Strip away the connections:**

``` r
remove_linkages(struc)
#> <glycan_structure[5]>
#> [1] Man(??-?)[Man(??-?)]Man(??-?)GlcNAc(??-?)GlcNAc(??-
#> [2] Gal(??-?)GalNAc(??-
#> [3] GlcNAc(??-?)[Gal(??-?)]GalNAc(??-
#> [4] Man(??-?)[Man(??-?)]Man(??-?)[Man(??-?)]Man(??-
#> [5] GlcNAc6Ac(??-?)Glc3Me(??-
#> # Unique structures: 5
```

**Remove the decorations:**

``` r
# Let's look at our decorated structure first
struc[5]
#> <glycan_structure[1]>
#> [1] GlcNAc6Ac(b1-4)Glc3Me(a1-
#> # Unique structures: 1

# Now remove the decorations (6Ac and 3Me)
remove_substituents(struc[5])
#> <glycan_structure[1]>
#> [1] GlcNAc(b1-4)Glc(a1-
#> # Unique structures: 1
```

**Convert monosaccharides to generic:**

``` r
convert_to_generic(struc)
#> <glycan_structure[5]>
#> [1] Hex(a1-3)[Hex(a1-6)]Hex(b1-4)HexNAc(b1-4)HexNAc(b1-
#> [2] Hex(b1-3)HexNAc(a1-
#> [3] Hex(b1-3)[HexNAc(b1-6)]HexNAc(a1-
#> [4] Hex(a1-3)[Hex(a1-6)]Hex(a1-3)[Hex(a1-6)]Hex(a1-
#> [5] HexNAc6Ac(b1-4)Hex3Me(a1-
#> # Unique structures: 5
```

**Reduce structure resolution level:**

``` r
reduce_structure_level(struc, to_level = "basic")
#> <glycan_structure[5]>
#> [1] Hex(??-?)[Hex(??-?)]Hex(??-?)HexNAc(??-?)HexNAc(??-
#> [2] Hex(??-?)HexNAc(??-
#> [3] HexNAc(??-?)[Hex(??-?)]HexNAc(??-
#> [4] Hex(??-?)[Hex(??-?)]Hex(??-?)[Hex(??-?)]Hex(??-
#> [5] HexNAc6Ac(??-?)Hex3Me(??-
#> # Unique structures: 5
# Same as remove_linkages() then convert_to_generic()
```

## Part 3: Conversions and Integrations

### From Structure to Composition

Ever wondered what’s actually in those complex structures? Easy:

``` r
comp <- as_glycan_composition(struc)
comp
#> <glycan_composition[5]>
#> [1] Man(3)GlcNAc(2)
#> [2] Gal(1)GalNAc(1)
#> [3] Gal(1)GlcNAc(1)GalNAc(1)
#> [4] Man(5)
#> [5] Glc(1)GlcNAc(1)Me(1)Ac(1)
```

### Back to Strings

Need to export your data or use it elsewhere?

``` r
# Get the original string representations
as.character(struc)
#> [1] "Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
#> [2] "Gal(b1-3)GalNAc(a1-"                                
#> [3] "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(a1-"                  
#> [4] "Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-"    
#> [5] "GlcNAc6Ac(b1-4)Glc3Me(a1-"
as.character(comp)
#> [1] "Man(3)GlcNAc(2)"           "Gal(1)GalNAc(1)"          
#> [3] "Gal(1)GlcNAc(1)GalNAc(1)"  "Man(5)"                   
#> [5] "Glc(1)GlcNAc(1)Me(1)Ac(1)"
```

### Playing Nice with the Tidyverse

`glyrepr` objects are first-class citizens in the tidyverse:

``` r
suppressPackageStartupMessages(library(tibble))
suppressPackageStartupMessages(library(dplyr))

df <- tibble(
  id = seq_along(struc),
  structures = struc,
  names = c("N-glycan core", "Core 1", "Core 2", "Branched Man", "Decorated")
)

df %>% 
  mutate(n_man = count_mono(structures, "Man")) %>%
  filter(n_man > 1)
#> # A tibble: 2 × 4
#>      id structures                                          names         n_man
#>   <int> <struct>                                            <chr>         <int>
#> 1     1 Man(a1-3)[Man(a1-6)]Man(b1-4)GlcNAc(b1-4)GlcNAc(b1- N-glycan core     3
#> 2     4 Man(a1-3)[Man(a1-6)]Man(a1-3)[Man(a1-6)]Man(a1-     Branched Man      5
```

## What’s Next?

Congratulations! You’ve just learned the fundamentals of glycan
representation in R. Here’s what you can explore next:

- 🔬 **Advanced analysis**: Check out the “Power User Guide: Efficient
  Glycan Manipulation” vignette for power-user features
- 🧬 **Motif searching**: Try the `glymotif` package for finding
  patterns in glycan structures  
- 📊 **Visualization**: Explore glycan visualization packages in the
  glycoverse

The glycoverse is your oyster! 🦪

## Session Information

``` r
sessionInfo()
#> R version 4.5.2 (2025-10-31)
#> Platform: x86_64-pc-linux-gnu
#> Running under: Ubuntu 24.04.3 LTS
#> 
#> Matrix products: default
#> BLAS:   /usr/lib/x86_64-linux-gnu/openblas-pthread/libblas.so.3 
#> LAPACK: /usr/lib/x86_64-linux-gnu/openblas-pthread/libopenblasp-r0.3.26.so;  LAPACK version 3.12.0
#> 
#> locale:
#>  [1] LC_CTYPE=C.UTF-8       LC_NUMERIC=C           LC_TIME=C.UTF-8       
#>  [4] LC_COLLATE=C.UTF-8     LC_MONETARY=C.UTF-8    LC_MESSAGES=C.UTF-8   
#>  [7] LC_PAPER=C.UTF-8       LC_NAME=C              LC_ADDRESS=C          
#> [10] LC_TELEPHONE=C         LC_MEASUREMENT=C.UTF-8 LC_IDENTIFICATION=C   
#> 
#> time zone: UTC
#> tzcode source: system (glibc)
#> 
#> attached base packages:
#> [1] stats     graphics  grDevices utils     datasets  methods   base     
#> 
#> other attached packages:
#> [1] dplyr_1.1.4        tibble_3.3.1       tictoc_1.2.1       glyrepr_0.9.0.9000
#> 
#> loaded via a namespace (and not attached):
#>  [1] jsonlite_2.0.0    compiler_4.5.2    tidyselect_1.2.1  stringr_1.6.0    
#>  [5] jquerylib_0.1.4   systemfonts_1.3.1 textshaping_1.0.4 yaml_2.3.12      
#>  [9] fastmap_1.2.0     R6_2.6.1          generics_0.1.4    igraph_2.2.1     
#> [13] knitr_1.51        backports_1.5.0   checkmate_2.3.3   rstackdeque_1.1.1
#> [17] desc_1.4.3        bslib_0.9.0       pillar_1.11.1     rlang_1.1.7      
#> [21] utf8_1.2.6        cachem_1.1.0      stringi_1.8.7     xfun_0.56        
#> [25] fs_1.6.6          sass_0.4.10       cli_3.6.5         pkgdown_2.2.0    
#> [29] magrittr_2.0.4    digest_0.6.39     lifecycle_1.0.5   vctrs_0.7.0      
#> [33] evaluate_1.0.5    glue_1.8.0        ragg_1.5.0        rmarkdown_2.30   
#> [37] purrr_1.2.1       tools_4.5.2       pkgconfig_2.0.3   htmltools_0.5.9
```
