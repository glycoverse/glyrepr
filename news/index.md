# Changelog

## glyrepr (development version)

## glyrepr 0.7.5

CRAN release: 2025-10-29

### Minor improvements and bug fixes

- IUPAC sequence generation is deterministic now for tie branches.
- Glycan structures are now truncated when printed to console as a
  column in a tibble.
- Glycan compositions are now colored when printed to console as a
  column in a tibble.

## glyrepr 0.7.4

CRAN release: 2025-09-23

### Minor improvements and bug fixes

- Update package title and description.
- Remove parallel examples in
  [`smap()`](https://glycoverse.github.io/glyrepr/reference/smap.md).

## glyrepr 0.7.3

### Minor improvements and bug fixes

- Fix some typos in the documentation.
- Add examples to some functions.
- Prepare for release on CRAN.

## glyrepr 0.7.2

### Minor improvements and bug fixes

- Fix the bug that
  [`structure_to_iupac()`](https://glycoverse.github.io/glyrepr/reference/structure_to_iupac.md)
  returns incorrect sequences with incorrect backbone or branch order.

## glyrepr 0.7.1

### Minor improvements and bug fixes

- Fix the bug that
  [`smap2()`](https://glycoverse.github.io/glyrepr/reference/smap2.md),
  [`spmap()`](https://glycoverse.github.io/glyrepr/reference/spmap.md),
  and related functions return unexpected results when the input `y` is
  a list.
- Improve the documentation of
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md),
  including the new behavior of vertex and edge order introduced in
  0.7.0.

## glyrepr 0.7.0

### Breaking changes

- `convert_mono_type()` is now replaced by
  [`convert_to_generic()`](https://glycoverse.github.io/glyrepr/reference/convert_to_generic.md).
  `convert_mono_type()` was created when three monosaccharide types
  existed: “concrete”, “generic”, and “simple”. When “simple” was
  removed, the old `convert_mono_type()` seems redundant, as the only
  valid conversion is from “concrete” to “generic” now. Therefore, we
  remove this function now and add a more straightforward
  [`convert_to_generic()`](https://glycoverse.github.io/glyrepr/reference/convert_to_generic.md).
- [`get_structure_graphs()`](https://glycoverse.github.io/glyrepr/reference/get_structure_graphs.md)
  is redesigned.
  - The `i` parameter is removed, as indexing can be done manually on
    the input `glyrepr_structure` vector or on the returned list easily.
  - Add a `return_list` parameter to control the return type. This
    parameter makes this function “type-stable”.
- [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  and
  [`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_structure.md)
  now reorder the underlying graphs to be in line with the IUPAC-style
  sequence. For example, the vertex order of
  “Gal(b1-3)\[GlcNAc(b1-6)\]GalNAc(b1-” is always 1. Gal, 2. GlcNAc, 3.
  GalNAc, and edges b1-3, b1-6, no matter what the original graphs are.
  Users can assign the indices of vertices and edges easily by printing
  the structure to console. This update makes `glymotif::match_motif()`
  more meaningful.

### New features

- `glyrepr_structure` has colors now in tibbles when printed to console.

### Minor improvements and bug fixes

- Fix the bug that monosaccharides with substituents are not colored
  when a
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  is printed in the console. For example, the “Neu5Ac” part in
  “Neu5Ac9Ac(a2-” was printed in black. Now it is printed in purple,
  while the “9Ac” part remains in black.
- Fix the bug that Neu5Ac and Neu5Gc with substituents at position 2, 3,
  or 4 could not be correctly parsed. Now, complex patterns like
  “Neu4Ac5Ac9Ac” can be properly handled, into a “Neu5Ac” monosaccharide
  with “4Ac,9Ac” as substituents.

## glyrepr 0.6.1

### Minor improvements and bug fixes

- [`n_glycan_core()`](https://glycoverse.github.io/glyrepr/reference/n_glycan_core.md)
  now has a “b1” reducing end anomer, not “?1”.
- Add validation to
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  to ensure no duplicated linkage positions. For example,
  “Gal(b1-3)\[Fuc(a1-3)\]GalNAc(b1-” is invalid now becuase both “Gal”
  and “Fuc” are linked to “GalNAc” at position 3.
- Add descriptions about ambiguous linkages and anomers in
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  documentation.
- [`remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.md)
  now also removes reducing end anomers.
- [`n_glycan_core()`](https://glycoverse.github.io/glyrepr/reference/n_glycan_core.md),
  [`o_glycan_core_1()`](https://glycoverse.github.io/glyrepr/reference/n_glycan_core.md),
  and
  [`o_glycan_core_2()`](https://glycoverse.github.io/glyrepr/reference/n_glycan_core.md)
  now have “??” anomers when `linkage = FALSE`.
- Fix a bug in
  [`smap_structure()`](https://glycoverse.github.io/glyrepr/reference/smap.md),
  [`smap2_structure()`](https://glycoverse.github.io/glyrepr/reference/smap2.md),
  [`spmap_structure()`](https://glycoverse.github.io/glyrepr/reference/spmap.md),
  and
  [`simap_structure()`](https://glycoverse.github.io/glyrepr/reference/simap.md)
  where modifying the structures can create identical structures, but
  the unique structures are not updated correctly. This automatically
  fixes a similar bug in
  [`remove_linkages()`](https://glycoverse.github.io/glyrepr/reference/remove_linkages.md).

## glyrepr 0.6.0

### Breaking changes

- Remove the “alditol” attribute from
  [`glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/glycan_structure.md)
  objects. This information is rarely used in glycomics and
  glycoproteomics data analysis. It is removed according to the razor
  principle.
- [`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_structure.md)
  now doesn’t allow the input IUPAC-condensed strings to omit the anomer
  information. Previously, something like “Glc(a1-3)GlcNAc” is valid.
  [`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/reference/as_glycan_structure.md)
  assumed that the core “GlcNAc” has a “?1-” anomer and added it
  automatically. The problem is that this behavior was not easily awared
  by users and might cause confusion. Again, less is more, so we remove
  it.

## glyrepr 0.5.2

### Minor improvements and bug fixes

- Fix some error message format errors.
- Update the abbreviated type name of `glyrepr_structure` from
  `structure` to `struct`.

## glyrepr 0.5.0

### Major changes

- Glycan structures now support multiple substituents on a single
  monosaccharide. Substituents are stored as comma-separated strings
  internally and concatenated in IUPAC format for display.

- Glycan compositions now support substituents. The `glycan_composition`
  class can now represent and count substituents alongside
  monosaccharides.
