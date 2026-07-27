# Enumerate Floating-Part Localizations

`enumerate_floating_localizations()` generates every conflict-free,
fully localized structure permitted by the candidate-parent domains in
`x`. Each variant records the complete assignment that produced it.

Candidate combinations are validated simultaneously, including linkages
with multiple possible acceptor positions such as `"a2-3/6"`. Variants
are canonicalized and then deduplicated by structure; when multiple
assignments produce the same canonical structure, the first assignment
in deterministic candidate order is retained.

`max_variants` is a conservative per-input safeguard. It limits the raw
Cartesian product before conflict filtering or canonical deduplication,
so the function may ask for a higher bound even when fewer variants
would ultimately remain.

Missing inputs and structures without floating parts each produce one
row with `variant_id = 1L`, the original structure, and an empty
assignment table. Empty structure vectors produce a zero-row result.

## Usage

``` r
enumerate_floating_localizations(x, max_variants = 256)
```

## Arguments

- x:

  A glycan structure vector.

- max_variants:

  A positive integer giving the maximum raw candidate combinations
  allowed for each input structure.

## Value

A tibble with columns:

- `input_id`: the integer position in `x`.

- `variant_id`: the sequential identifier after canonical deduplication.

- `structure`: a `glyrepr_structure` vector column containing fully
  localized variants.

- `assignments`: a list-column of tibbles with `glycan_id`, `part_id`,
  and `parent_node`. Here `glycan_id` equals `input_id`.

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|1,2}Gal(b1-3)GalNAc(a1-"
)
enumerate_floating_localizations(glycan)
#> # A tibble: 2 × 4
#>   input_id variant_id structure                         assignments     
#>      <int>      <int> <struct>                          <list>          
#> 1        1          1 Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-   <tibble [1 × 3]>
#> 2        1          2 Gal(b1-3)[Neu5Ac(a2-6)]GalNAc(a1- <tibble [1 × 3]>
```
