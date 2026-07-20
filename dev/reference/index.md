# Package index

## Glycan Composition and Structure Creation

- [`glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_composition.md)
  [`is_glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_composition.md)
  : Create a Glycan Composition
- [`as_glycan_composition()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_composition.md)
  : Convert to Glycan Composition
- [`glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
  [`is_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/glycan_structure.md)
  : Create a Glycan Structure Vector
- [`as_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/as_glycan_structure.md)
  : Convert to Glycan Structure Vector
- [`structure_nodes()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
  [`structure_edges()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
  [`structure_from_tibbles()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
  : Convert Glycan Structures to Graph Tables

## Glycan Composition and Structure Inspection

- [`count_mono()`](https://glycoverse.github.io/glyrepr/dev/reference/count_mono.md)
  : Get the Number of Monosaccharides
- [`get_anomer()`](https://glycoverse.github.io/glyrepr/dev/reference/get_anomer.md)
  : Get the Anomeric information
- [`get_mono_type()`](https://glycoverse.github.io/glyrepr/dev/reference/get_mono_type.md)
  : Get Monosaccharide Types
- [`get_structure_graphs()`](https://glycoverse.github.io/glyrepr/dev/reference/get_structure_graphs.md)
  : Access Individual Glycan Structures
- [`has_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/has_linkages.md)
  : Determine if a Glycan Structure has Linkages
- [`structure_nodes()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
  [`structure_edges()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
  [`structure_from_tibbles()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_tables.md)
  : Convert Glycan Structures to Graph Tables
- [`structure_to_iupac()`](https://glycoverse.github.io/glyrepr/dev/reference/structure_to_iupac.md)
  : Convert Glycan Structure to IUPAC-like Sequence
- [`get_structure_level()`](https://glycoverse.github.io/glyrepr/dev/reference/get_structure_level.md)
  : Get the Structure Resolution Levels

## Low-Level Glycan Graph Infrastructure

Developer-facing functions for validating and canonicalizing glycan
graphs, generating IUPAC-condensed keys, and assembling trusted
structure vectors.

- [`validate_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph.md)
  : Validate a Glycan Graph
- [`canonicalize_glycan_graph()`](https://glycoverse.github.io/glyrepr/dev/reference/canonicalize_glycan_graph.md)
  : Canonicalize a Glycan Graph
- [`validate_glycan_graph_vector()`](https://glycoverse.github.io/glyrepr/dev/reference/validate_glycan_graph_vector.md)
  : Validate Compatibility Across Glycan Graphs
- [`graph_to_iupac()`](https://glycoverse.github.io/glyrepr/dev/reference/graph_to_iupac.md)
  : Generate IUPAC-Condensed from a Glycan Graph
- [`new_glycan_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/new_glycan_structure.md)
  : Construct a Glycan Structure Vector from Trusted Data

## Glycan Composition and Structure Manipulation

- [`convert_to_generic()`](https://glycoverse.github.io/glyrepr/dev/reference/convert_to_generic.md)
  : Convert Monosaccharides to Generic Type
- [`fill_anomer_pos()`](https://glycoverse.github.io/glyrepr/dev/reference/fill_anomer_pos.md)
  : Fill Anomer Positions
- [`remove_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/remove_linkages.md)
  : Remove All Linkages from a Glycan
- [`remove_substituents()`](https://glycoverse.github.io/glyrepr/dev/reference/remove_substituents.md)
  : Remove All Substituents from a Glycan
- [`reduce_structure_level()`](https://glycoverse.github.io/glyrepr/dev/reference/reduce_structure_level.md)
  : Reduce a Glycan Structure to a Lower Resolution Level
- [`simap()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  [`simap_vec()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  [`simap_lgl()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  [`simap_int()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  [`simap_dbl()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  [`simap_chr()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  [`simap_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/simap.md)
  : Map Functions Over Glycan Structure Vectors with Indices
- [`smap()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  [`smap_vec()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  [`smap_lgl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  [`smap_int()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  [`smap_dbl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  [`smap_chr()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  [`smap_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/smap.md)
  : Map Functions Over Glycan Structure Vectors
- [`smap2()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  [`smap2_vec()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  [`smap2_lgl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  [`smap2_int()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  [`smap2_dbl()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  [`smap2_chr()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  [`smap2_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/smap2.md)
  : Map Functions Over Two Glycan Structure Vectors
- [`ssome()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_predicates.md)
  [`severy()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_predicates.md)
  [`snone()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_predicates.md)
  : Test Predicates on Glycan Structure Vectors
- [`smap_unique()`](https://glycoverse.github.io/glyrepr/dev/reference/smap_unique.md)
  : Apply Function to Unique Structures Only
- [`spmap()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  [`spmap_vec()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  [`spmap_lgl()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  [`spmap_int()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  [`spmap_dbl()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  [`spmap_chr()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  [`spmap_structure()`](https://glycoverse.github.io/glyrepr/dev/reference/spmap.md)
  : Map Functions Over Glycan Structure Vectors and Multiple Arguments

## Example Glycan Structures

- [`n_glycan_core()`](https://glycoverse.github.io/glyrepr/dev/reference/n_glycan_core.md)
  [`o_glycan_core_1()`](https://glycoverse.github.io/glyrepr/dev/reference/n_glycan_core.md)
  [`o_glycan_core_2()`](https://glycoverse.github.io/glyrepr/dev/reference/n_glycan_core.md)
  : Example Glycan Structures

## Helper Functions

- [`available_monosaccharides()`](https://glycoverse.github.io/glyrepr/dev/reference/available_monosaccharides.md)
  : Get Available Monosaacharides
- [`available_substituents()`](https://glycoverse.github.io/glyrepr/dev/reference/available_substituents.md)
  : Available Substituents
- [`infer_anomer_pos()`](https://glycoverse.github.io/glyrepr/dev/reference/infer_anomer_pos.md)
  [`get_anomer_pos()`](https://glycoverse.github.io/glyrepr/dev/reference/infer_anomer_pos.md)
  : Infer Anomer Positions
- [`is_known_monosaccharide()`](https://glycoverse.github.io/glyrepr/dev/reference/is_known_monosaccharide.md)
  : Check if a Monosaccharide is Known
- [`possible_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/possible_linkages.md)
  : Generate Possible Linkages
- [`valid_linkages()`](https://glycoverse.github.io/glyrepr/dev/reference/valid_linkages.md)
  : Check if Linkages are Valid
