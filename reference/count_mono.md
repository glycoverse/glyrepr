# Get the Number of Monosaccharides

When `mono` is:

- `NULL` (default), returns the total number of monosaccharides and
  substituents.

- A string, returns the number of the specified monosaccharide or
  substituent.

## Usage

``` r
count_mono(x, mono = NULL, include_subs = FALSE)

# S3 method for class 'glyrepr_composition'
count_mono(x, mono = NULL, include_subs = FALSE)

# S3 method for class 'glyrepr_structure'
count_mono(x, mono = NULL, include_subs = FALSE)
```

## Arguments

- x:

  A glycan composition (`glyrepr_composition`) or a glycan structure
  (`glyrepr_structure`) vector.

- mono:

  The monosaccharide or substituent to count. A character scalar. If
  `NULL` (default), return the total number of monosaccharides.

- include_subs:

  Whether to include substituents when `mono` is `NULL`. Default is
  `FALSE`.

## Value

A numeric vector of the same length as `x`.

## Details

When `mono` is "generic" (e.g. "Hex", "HexNAc"), it counts all
"concrete" monosaccharides that match. For example, "Hex" will count all
Glc, Man, Gal, etc. When `mono` is "concrete" (e.g. "Gal", "GalNAc"), NA
is returned when the composition is "generic".

## Examples

``` r
comp <- glycan_composition(c(Gal = 1, Man = 1, GalNAc = 1))
count_mono(comp, "Hex")
#> [1] 2
count_mono(comp, "Gal")
#> [1] 1

struct <- as_glycan_structure("Gal(b1-3)GlcNAc(b1-4)Glc(a1-")
count_mono(struct, "Gal")
#> [1] 1

# Total number of monosaccharides
count_mono(comp)
#> [1] 3
```
