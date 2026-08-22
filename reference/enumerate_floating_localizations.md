# Enumerate Floating Localizations

`enumerate_floating_localizations()` generates every conflict-free,
fully localized structure permitted by the candidate-parent domains in
`x`. Each variant records the complete assignment that produced it.

Candidate combinations for floating parts and substituents are validated
simultaneously, including linkages or substituents with multiple
possible carbon positions. Floating-component dependencies must be
acyclic and ultimately connect to the main tree. Variants are
canonicalized and, by default, deduplicated by structure. When multiple
assignments produce the same canonical structure, the first assignment
in deterministic candidate order is retained. Set `deduplicate = FALSE`
to retain every valid assignment and its original-node provenance,
including assignments that produce identical canonical structures.

`max_variants` is a conservative per-input safeguard. It limits the raw
Cartesian product before conflict filtering or canonical deduplication,
so the function may ask for a higher bound even when fewer variants
would ultimately remain.

Missing inputs and structures without floating metadata each produce one
row with `variant_id = 1L`, the original structure, and an empty
assignment table. Empty structure vectors produce a zero-row result.

## Usage

``` r
enumerate_floating_localizations(x, max_variants = 256, deduplicate = TRUE)
```

## Arguments

- x:

  A glycan structure vector.

- max_variants:

  A positive integer giving the maximum raw candidate combinations
  allowed for each input structure.

- deduplicate:

  A logical value. If `TRUE` (default), retain only the first assignment
  for each canonical structure. If `FALSE`, retain every conflict-free
  assignment.

## Value

A tibble with columns:

- `input_id`: the integer position in `x`.

- `variant_id`: the sequential identifier after optional canonical
  deduplication.

- `structure`: a `glyrepr_structure` vector column containing fully
  localized variants.

- `assignments`: a list-column of tibbles with `glycan_id`, `part_id`,
  `parent_node`, and `substituent_id`. Exactly one of `part_id` and
  `substituent_id` is non-missing in each row. Here `glycan_id` equals
  `input_id`.

## Examples

``` r
glycan <- as_glycan_structure(
  "{Neu5Ac(a2-6)|2,3}Gal(b1-3)GalNAc(a1-"
)
enumerate_floating_localizations(glycan)
#> # A tibble: 2 × 4
#>   input_id variant_id structure                         assignments     
#>      <int>      <int> <struct>                          <list>          
#> 1        1          1 Neu5Ac(a2-6)Gal(b1-3)GalNAc(a1-   <tibble [1 × 4]>
#> 2        1          2 Gal(b1-3)[Neu5Ac(a2-6)]GalNAc(a1- <tibble [1 × 4]>
```
