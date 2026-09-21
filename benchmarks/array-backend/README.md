# Parser backend integration

The implementation adds two public R interfaces. Other packages do not need to
link to private C++ symbols or depend on igraph's internal representation.

## New parsers: residue and edge arrays

`structure_from_arrays(records, on_failure = "error")` accepts a list of records.
Each record describes one glycan; `NULL` represents a missing glycan. Required
fields are `mono`, `sub`, `edges`, `linkage`, and `anomer`. Optional fields are
`alditol`, `floating_parts`, and `floating_substituents`.

Node IDs are one-based array positions. `edges` contains alternating parent and
child IDs, not a matrix. Floating metadata uses those same IDs, including
cross-component candidates. See `?structure_from_arrays` for the full schema.

```r
records <- list(example = list(
  mono = c("Glc", "Gal"), sub = c("", ""),
  edges = c(1L, 2L), linkage = "b1-4", anomer = "?1"
))
x <- glyrepr::structure_from_arrays(records)
```

The R boundary checks record types and chemical tokens. C++ checks topology,
component membership, candidate feasibility, and linkage conflicts and computes
canonical order and IUPAC keys. Successful records only materialize an igraph
after canonicalization, once per distinct canonical key. Floating symmetry uses
a public igraph BLISS callback. Source-array order need not be canonical.

Both the IUPAC parser and array importer use `src/compact-core.h`. Its internal
representation is not a public ABI. The native 2048-node guard selects the R graph
reference path; it is not a public size limit. The reference implementation's
own recursion and resource limits still apply.

## Existing graph parsers: combined canonicalization

`canonicalize_glycan_graphs(graphs)` returns aligned `graphs`, `iupac`, `status`,
and `reason` fields. It shares ordering/key traversal for ordinary trees and
uses the established R graph path for floating metadata. It does not yet move
arbitrary igraph canonicalization into C++.

It deliberately retains one graph per input, including duplicate identities with
different source attributes. The caller can then deduplicate explicitly when
constructing a structure vector:

```r
batch <- glyrepr::canonicalize_glycan_graphs(graphs, on_failure = "na")
keep <- batch$status == "ok" & !duplicated(batch$iupac)
x <- glyrepr::new_glycan_structure(
  batch$iupac,
  stats::setNames(batch$graphs[keep], batch$iupac[keep])
)
```

Only pass `validate = FALSE` for individually validated graphs. This preserves
the existing trusted-input contract. The graph interface retains arbitrary
vertex, edge, and graph attributes through the established graph permutations.

Strict failures from both new interfaces have class
`glyrepr_error_structure_failure`, with `position`, `input_name`, and `reason`.
Recovery warnings have class `glyrepr_warning_structure_failure`, with
`positions`, `reasons`, and `input_names`. Recovery does not discard successful
records or turn missing input into a failure.

## Validation and reproduction

Run from the package root:

```sh
Rscript benchmarks/array-backend/verify.R
R CMD INSTALL --preclean --library=/path/to/release/library .
Rscript benchmarks/array-backend/benchmark.R /path/to/release/library
Rscript benchmarks/array-backend/check-reference-roundtrip.R
```

`verify.R` uses the frozen corpus and extended inputs already tracked under
`benchmarks/compact-cpp/results-v2`. It parses with the R reference parser,
reverses source node IDs (including floating metadata), and compares the new
constructor against the R reference's complete canonical graph attributes and
key. The final audit covers 11,425 distinct strings: 10,640 valid records all
have native status `ok` and exact canonical-key / semantic-graph parity; 785
reference-invalid strings are recorded separately, not counted as valid parity.
No valid records needed graph fallback in this audit.

The benchmark compares equivalent complete structure-vector outputs and checks
key equality before timing. Graph split/fused timings both include validation.
`native_arrays_only` excludes record validation and graph materialization and
must not be interpreted as complete constructor speed. Each workload has three
alternating-order measurement rounds after warmup.

## Existing reference round-trip limitation

During the audit, converting *already canonicalized reference graphs* back to
array input exposed two pre-existing floating normalization cases. Resolving a
singleton floating attachment can leave other explicit candidate sets pointing
to now-occupied carbon positions. The R validator then rejects its own canonical
output. The corresponding raw parser records are valid and the new array backend
produces the same output as the R reference.

The inputs, reference graphs, and reproduction output are preserved under
`known-reference/`. This change does not modify that existing candidate-domain
normalization policy. These cases remain an explicit limitation, not a claim of
successful canonical-graph revalidation.

The glyparse source has not been migrated in this change. It can adopt the graph
batch interface first, then emit array records directly from individual parsers.
`auto_parse()` grouping and format-specific residue caches remain glyparse work.
