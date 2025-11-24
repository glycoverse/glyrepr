# glyrepr (development version)

# glyrepr 0.9.0

## New features

* Update the monosaccharides table, ensuring all monosaccharides having generic names:
  * Add generic monosaccharides to Ara, Lyx, Xyl, Rib, Api, Neu, Kdn, Pse, Leg, Aci, 4eLeg, Bac, LDmanHep, DDmanHep, Kdo, Dha, MurNAc, MurNGc, Mur, Fru, Tag, Sor, Psi.
  * Rename "Pent" to "Pen".
  * Delete "Sia" from the table.
* `count_mono()` now supports counting substituents, with a new argument `include_subs`.
* Add `reduce_structure_level()` to reduce a glycan structure to a lower resolution level.

## Minor improvements and bug fixes

* Fix the bug that `convert_to_generic()` fails with glycan compositions containing substituents.
* Fix the bug that `has_linkages()` didn't consider reducing end anomers, which influenced the results of `get_structure_level()`.

# glyrepr 0.8.0

## Breaking changes

* Remove `normalize_substituents()`. This function does not help outside of `glyrepr` so we make it an internal function.
* `glycan_composition()` cannot accept empty integer vectors now. Therefore, `glycan_composition(integer(0))` is not valid anymore.
* `glycan_composition()` now checks input types more strictly.
* `as_glycan_composition()` handles NA values and empty strings ("") more strictly. Now an error will be raised if NA values or empty strings are passed instead of dropping them silently. This update makes `as_glycan_composition()` size-stable and consistent with `as_glycan_structure()`.

## New features

* Add `get_structure_level()` to get the structure resolution levels of a glycan structure vector.
* `as_glycan_composition()` now supports simple composition strings like "H5N2", "H5N4F1S2", "H5N4A1G1", etc.
* `count_mono()` now returns total number of monosaccharides when `mono` is `NULL`.
* `has_linkages()` now has a `strict` parameter to control the strictness of the check.

## Minor improvements and bug fixes

* `glycan_composition()` now supports `!!!`.
* Add more examples about character strings in the documentation of `as_glycan_composition()`.
* Add a section about structure resolution levels in the getting started vignette.

# glyrepr 0.7.5

## Minor improvements and bug fixes

* IUPAC sequence generation is deterministic now for tie branches.
* Glycan structures are now truncated when printed to console as a column in a tibble.
* Glycan compositions are now colored when printed to console as a column in a tibble.

# glyrepr 0.7.4

## Minor improvements and bug fixes

* Update package title and description.
* Remove parallel examples in `smap()`.

# glyrepr 0.7.3

## Minor improvements and bug fixes

* Fix some typos in the documentation.
* Add examples to some functions.
* Prepare for release on CRAN.

# glyrepr 0.7.2

## Minor improvements and bug fixes

* Fix the bug that `structure_to_iupac()` returns incorrect sequences with incorrect backbone or branch order.

# glyrepr 0.7.1

## Minor improvements and bug fixes

* Fix the bug that `smap2()`, `spmap()`, and related functions return unexpected results when the input `y` is a list.
* Improve the documentation of `glycan_structure()`, including the new behavior of vertex and edge order introduced in 0.7.0.

# glyrepr 0.7.0

## Breaking changes

* `convert_mono_type()` is now replaced by `convert_to_generic()`. `convert_mono_type()` was created when three monosaccharide types existed: "concrete", "generic", and "simple". When "simple" was removed, the old `convert_mono_type()` seems redundant, as the only valid conversion is from "concrete" to "generic" now. Therefore, we remove this function now and add a more straightforward `convert_to_generic()`.
* `get_structure_graphs()` is redesigned.
  - The `i` parameter is removed, as indexing can be done manually on the input `glyrepr_structure` vector or on the returned list easily.
  - Add a `return_list` parameter to control the return type. This parameter makes this function "type-stable".
* `glycan_structure()` and `as_glycan_structure()` now reorder the underlying graphs to be in line with the IUPAC-style sequence. For example, the vertex order of "Gal(b1-3)[GlcNAc(b1-6)]GalNAc(b1-" is always 1. Gal, 2. GlcNAc, 3. GalNAc, and edges b1-3, b1-6, no matter what the original graphs are. Users can assign the indices of vertices and edges easily by printing the structure to console. This update makes `glymotif::match_motif()` more meaningful. 

## New features

* `glyrepr_structure` has colors now in tibbles when printed to console.

## Minor improvements and bug fixes

* Fix the bug that monosaccharides with substituents are not colored when a `glycan_structure()` is printed in the console. For example, the "Neu5Ac" part in "Neu5Ac9Ac(a2-" was printed in black. Now it is printed in purple, while the "9Ac" part remains in black.
* Fix the bug that Neu5Ac and Neu5Gc with substituents at position 2, 3, or 4 could not be correctly parsed. Now, complex patterns like "Neu4Ac5Ac9Ac" can be properly handled, into a "Neu5Ac" monosaccharide with "4Ac,9Ac" as substituents.

# glyrepr 0.6.1

## Minor improvements and bug fixes

* `n_glycan_core()` now has a "b1" reducing end anomer, not "?1".
* Add validation to `glycan_structure()` to ensure no duplicated linkage positions. For example, "Gal(b1-3)[Fuc(a1-3)]GalNAc(b1-" is invalid now becuase both "Gal" and "Fuc" are linked to "GalNAc" at position 3.
* Add descriptions about ambiguous linkages and anomers in `glycan_structure()` documentation.
* `remove_linkages()` now also removes reducing end anomers.
* `n_glycan_core()`, `o_glycan_core_1()`, and `o_glycan_core_2()` now have "??" anomers when `linkage = FALSE`.
* Fix a bug in `smap_structure()`, `smap2_structure()`, `spmap_structure()`, and `simap_structure()` where modifying the structures can create identical structures, but the unique structures are not updated correctly. This automatically fixes a similar bug in `remove_linkages()`.

# glyrepr 0.6.0

## Breaking changes

* Remove the "alditol" attribute from `glycan_structure()` objects. This information is rarely used in glycomics and glycoproteomics data analysis. It is removed according to the razor principle.
* `as_glycan_structure()` now doesn't allow the input IUPAC-condensed strings to omit the anomer information. Previously, something like "Glc(a1-3)GlcNAc" is valid. `as_glycan_structure()` assumed that the core "GlcNAc" has a "?1-" anomer and added it automatically. The problem is that this behavior was not easily awared by users and might cause confusion. Again, less is more, so we remove it.

# glyrepr 0.5.2

## Minor improvements and bug fixes

* Fix some error message format errors.
* Update the abbreviated type name of `glyrepr_structure` from `structure` to `struct`.

# glyrepr 0.5.0

## Major changes

* Glycan structures now support multiple substituents on a single monosaccharide.
  Substituents are stored as comma-separated strings internally and concatenated
  in IUPAC format for display.

* Glycan compositions now support substituents. The `glycan_composition` class
  can now represent and count substituents alongside monosaccharides.
