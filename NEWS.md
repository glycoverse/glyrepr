# glyrepr 1.0.0

## Breaking changes

* `get_mono_type()` now returns one value per structure or composition instead of one scalar for the whole vector and may return `"mixed"`; callers that assume a scalar result must handle aligned per-element values. Generic and concrete residues can now coexist within structures, compositions, and their vectors. (#90)
* Structure levels now depend only on linkage and anomer information. `get_structure_level()` returns one `"intact"`, `"partial"`, or `"topological"` value per structure; the former `"basic"` level and `reduce_structure_level()` have been removed. Use `remove_linkages()` to remove linkage information and `convert_to_generic()` to convert residue identities. (#90)

## New features

* Generic and concrete residues can now coexist within structures, compositions, and their vectors; use aligned `get_mono_type()` results to inspect each element. (#90)
* Glycan structures now support alditols through reducing-end `-ol` IUPAC syntax, such as `Gal(b1-4)GlcNAc-ol(a1-`. New `get_alditol()` inspects the graph-level status, and `structure_from_tibbles()` gains an `alditols` argument for lossless table round-trips; legacy graphs without the attribute remain non-alditols. (#88)
* Concrete monosaccharides now support explicit unusual absolute configurations using a leading `D-` or `L-`, such as `D-Fuc`, `L-Gul`, and `D-Fucf`; unprefixed names retain their natural configurations. (#85, #86)
* Concrete monosaccharides now support explicit furanose forms such as `Galf` and `GlcfNAc`; generic conversion remains unchanged, so these become `Hex` and `HexNAc`, respectively. (#82)
* `as_glycan_structure()` and `glycan_structure()` now support floating glycan substructures with optional candidate-parent indices using `{<floating>}` and `{<floating>|<parents>}` syntax. Floating structures retain component-node and candidate-parent metadata through validation, canonicalization, transformations, and graph-table round-trips. (#80)
* `as_glycan_structure()` and `glycan_structure()` now support floating substituents with unresolved parent residues using `{<substituent>}` and `{<substituent>|<parents>}` syntax, including unknown carbon positions such as `{?S}`. (#87)
* Floating candidate-parent indices follow complete IUPAC-condensed residue order. Floating parts may attach to other floating components or the main tree, floating substituents may target any residue node, and localization retains only conflict-free acyclic assignments that connect every component to the main tree. (#89)
* New `enumerate_floating_graph_localizations()` and `enumerate_floating_localizations()` return conflict-free floating-part and floating-substituent assignments with provenance, while `localize_floating_parts()` attaches selected floating parts; graph-level results retain original vertex IDs, and `deduplicate = FALSE` retains assignments that canonicalize to the same structure. (#80)
* `print()` gains an `n` argument for `glyrepr_structure` and `glyrepr_composition` vectors. (#77)
* New `structure_candidate_edges()`, `structure_component_membership()`, and `structure_floating_candidates()` helpers expose floating-part membership and virtual edges plus floating-part and floating-substituent candidate parents for inspection, drawing, and constraint-aware graph operations. (#80)
* Structure inspection, transformation, graph-table, composition, and IUPAC APIs now accept individual glycan `igraph` objects directly while preserving existing `glyrepr_structure` vector behavior and vertex IDs. (#81)

## Minor improvements and fixes

* `as_glycan_structure()` and `glycan_structure()` now support ambiguous substituent positions such as `Gal4/6S(a1-` and preserve them during graph and IUPAC conversion. (#84)
* `as_glycan_structure()` now accepts omitted reducing-end annotations by inferring the anomer position, normalizes `(?-?)` to `(??-?)`, and collapses linkage-position choices containing `?` to a single unknown position. (#83)
* `get_structure_level()` treats floating candidate-parent ambiguity independently from linkage resolution; fully specified floating structures can be intact. (#80)
* IUPAC-condensed sequence generation from glycan graphs is now faster. (#78)
* `validate_glycan_graph()` and `as_glycan_structure()` are faster through bulk linkage-position validation. (#79)

# glyrepr 0.14.0

## New features

* The low-level `validate_glycan_graph()`, `canonicalize_glycan_graph()`, `validate_glycan_graph_vector()`, `graph_to_iupac()`, and `new_glycan_structure()` APIs support name-preserving construction from trusted glycan graphs. (#75)

## Minor improvements and bug fixes

* Glycan structures now allow multiple substituents with unknown positions. (#67)
* `as_glycan_structure()` gains `on_failure = "na"` to preserve valid elements, replace element-local failures with `NA`, and report one aggregated warning. (#73)
* `convert_to_generic()` now converts concrete monosaccharides in mixed concrete/generic character vectors. (#68)
* `fill_anomer_pos()` now accepts glycan structures with generic monosaccharides. (#70)

# glyrepr 0.13.0

## Breaking changes

* Remove the `.parallel` argument from `smap()`, `smap2()`, `spmap()`, and their variants, and drop the unused `furrr` and `future` dependencies. (#62)

## New features

* Add `structure_nodes()`, `structure_edges()`, and `structure_from_tibbles()` for converting glycan structures to and from graph-table tibbles. (#60)
* Rename `get_anomer_pos()` to `infer_anomer_pos()` to better describe that the position is inferred from a monosaccharide name; `get_anomer_pos()` remains available as a backward-compatible alias. (#61)
* `get_anomer_pos()` now accepts generic monosaccharide names such as `"Hex"`. (#57)
* Add `NGc` and `Gc` to the supported substituent list for N-glycolyl and glycolyl substituents. (#54)
* Add `as.list()` support for `glyrepr_composition` and `glyrepr_structure` vectors. (#52)

## Minor improvements and bug fixes

* `glyrepr_structure` objects no longer return `TRUE` from `is.character()`. Use `as.character()` for explicit IUPAC-condensed string conversion. (#59)
* Fix `as_glycan_structure()` so character vectors containing `NA` reject mixed concrete and generic structures consistently. (#58)
* Fix `as_glycan_structure()` parsing for monosaccharide names that start with digits, such as `6dGul` and `4eLeg`. (#55)
* Fix `as_glycan_structure()` parsing for substituent names such as `Pyr`, `PC`, `PPEtn`, and `PEtn` that share the `P` prefix. (#53)
* Ambiguous linkages like `a2-3/6` are also regarded as unknown for `has_linkage()` and `get_structure_level()`. (#51)

# glyrepr 0.12.1

## Minor improvements and bug fixes

* Performance optimization for glycan structure vector creation. (#46)

# glyrepr 0.12.0

## New features

* Add `get_anomer_pos()` helper to get the anomer position of a monosaccharide.
* Add `fill_anomer_pos()` function to fill the anomer position of a glycan structure with missing anomer information.

# glyrepr 0.11.0

## Breaking changes

* `get_structure_level()` now returns one character scalar for a `glyrepr_structure` vector instead of one value per element. The vector-wide level is "intact", "partial", "topological", or "basic" according to the combined residue and linkage detail of the non-missing structures in the vector (#42).

## New features

* `as_glycan_composition()` now supports parsing "E" and "L" in the input composition strings as "NeuAc". For example, `as_glycan_composition("H5N4F1L1E1")` is now correctly parsed as `Hex(5)HexNAc(4)Fuc(1)NeuAc(2)`, with a warning about dropping the sialic acid linkage information (#41).

## Minor improvements and bug fixes

* Fix the bug that `glycan_composition()` and `as_glycan_composition()` cannot handle duplications in the input. For example, `as_glycan_composition("Hex(2)Hex(1)HexNAc(2)")` is correctly regared as `Hex(3)HexNAc(2)` now (#40).
* `as_glycan_structure(NA_character_)` now creates a missing structure instead of erroring.
* `get_structure_level()` now ignores missing structures when determining the vector-wide level, and returns `NA_character_` for empty or all-missing structure vectors.
* `reduce_structure_level()` preserves missing structures in output.
* `simap()` and ``simap_structure()` now skip missing structures like the other smap variants.
* `get_mono_type.glyrepr_composition()` now ignores missing composition elements and returns `NA_character_` for all-NA composition vectors.
* Rewrite the Getting Started vignette for better readability and adopt a calmer tone for all vignettes.

# glyrepr 0.10.1

## Minor improvements and bug fixes

* Use `dplyr::recode_values()` to replace deprecated `dplyr::case_match()` to prevent warnings from `dplyr`.

# glyrepr 0.10.0

We have redesigned the internal implementation of `glyrepr_composition` and `glyrepr_structure`. This brought native support for names to `glyrepr_structure`, and NA values to both `glyrepr_structure` and `glyrepr_composition`.

## New features

* `smap()`, `smap2()`, `spmap()`, `simap()` and their variants now preserve
  names from input `glyrepr_structure` vectors in their output.
* `glyrepr_structure` now formally supports names. All operations on a named `glyrepr_structure` vectors preserve the names.
* NA values are supported for `glyrepr_structure` and `glyrepr_composition`. Any operation on a `glyrepr_structure` or `glyrepr_composition` vector with NA values behave intuitively. `is.na()` now works for these two classes.
* `glycan_composition()` now accepts another `glyrepr_composition` vector as input, returning it as-is.

## Breaking changes

* `glyrepr_composition` and `glyrepr_structure` now enforce the same monosaccharide type ("concrete" or "generic") within a vector. Mixed types are not allowed anymore. This invariant is enforced both when creating new vectors and when combining existing vectors.
* `glycan_structure()` now does not support multiple `glyrepr_structure` vectors as input anymore. For example, `glycan_structure(o_glycan_core_1(), o_glycan_core_2())` is not valid anymore. Please use `c(o_glycan_core_1(), o_glycan_core_2())` instead.
* `get_mono_type()` now returns a character scalar instead of a character vector for `glyrepr_structure` and `glyrepr_composition`.

## Minor improvements and bug fixes

* Subsetting `glyrepr_structure` with `integer(0)` and `NULL` correctly removes all underlying graphs.
* `[[<-` is forbidden on `glyrepr_structure` vectors. Previously, the operation could be performed silently, but resulted in an invalid object.

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
