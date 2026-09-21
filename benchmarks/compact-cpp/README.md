# Compact C++ parsing prototype

This experiment leaves package code and dependencies unchanged. It tests whether
parsing and canonicalizing **before** materializing an igraph removes the graph
allocation/permutation bottleneck. It is not a production replacement.

## Design

`compact.cpp` scans IUPAC right to left into residue, linkage, parent-index, and
child-index arrays. Parents already exist when children are parsed. It validates
syntax, residue membership, linkage syntax, tree construction, and duplicate
concrete acceptor positions on these arrays. It computes subtree depth/signature,
selects the backbone, then emits canonical IUPAC plus final vertex/edge arrays in
one traversal. The R wrapper constructs the topology once with
`igraph::make_graph()`, then attaches attributes in batches. It never permutes or
rebuilds the graph, and C++ never calls igraph or reads its internal object layout.
Canonical duplicates share one stored graph; input names, duplicates and NA are
preserved.

The bounded fast path accepts known, unmodified residues (including concrete,
generic and unusual configurations), branches, ambiguous linkages, reducing-end
annotations, inferred donor position, and alditols. Modified residues and floating
metadata go to the existing public constructor. This preserves their validation
and canonicalization; it is **not evidence of C++ support for those features**.
Unknown residues and errors also return to the public constructor to preserve
its diagnostics. The prototype implements the strict character-input constructor;
it does not implement the `on_failure = "na"` interface or graph inputs.

R's active collation matters for tied branch signatures. Under `LC_COLLATE=C`,
the kernel sorts strings directly. Under other locales, it calls `base::order()`
on small signature vectors; all graph/tree work remains on arrays. Thus the
normal-locale prototype is not entirely free of R callbacks. Residue vocabulary
and donor positions are initialized from package tables, never from the benchmark
corpus. Strings are parsed afresh on every call (deduplicated within the call).
The 2,048-node guard routes oversized trees to the reference implementation.

## Reproduce

From the package root, with package development dependencies, Rcpp and glydb:

```sh
Rscript benchmarks/compact-cpp/verify.R
Rscript benchmarks/compact-cpp/audit-syntax.R
Rscript benchmarks/compact-cpp/benchmark.R 1
Rscript benchmarks/compact-cpp/benchmark.R 2
Rscript benchmarks/compact-cpp/benchmark.R 3
Rscript benchmarks/compact-cpp/audit-graph-calls.R
Rscript benchmarks/compact-cpp/summarize.R
```

Run workers sequentially to avoid benchmark contention. Each starts in a fresh R
process; compilation/source time is recorded separately and excluded from steady
state timings. Each worker runs two paired repetitions with alternating order,
explicit garbage collection before measurements, and three-element warmups.
Seed 20260921 freezes workloads. Full constructor timings include parsing,
validation, canonicalization, graph creation, deduplication, vector construction,
and any fallback. This is compared with the current checkout's public R API.

Isolated `parse_validate_canonical_pairs` returns equivalent graph/IUPAC pairs in
both implementations without vector assembly/deduplication. `arrays_only` and
`graph_materialization` report absolute costs for different stages; they are not
compared as if their outputs were equivalent. CPU contention and GC variation
remain possible; retain ranges and paired ratios rather than extrapolating one
measurement to the whole package.

`results/` contains corpus/workload snapshots, row-level fast-path coverage, exact
parity outcomes, raw repeated timings, paired speedups, session information and
sampled profiles. Verification compares canonical strings, exact ordered edge
endpoints, all vertex/edge/graph attributes, graph dictionary keys, names, missing
values and canonical deduplication. It includes all glydb entries, every known
residue and its alditol, 500 seeded noncanonical random trees, malformed inputs,
multiple collation settings, and literals extracted from existing parser tests.

See `RESULTS.md` for measured outcomes and remaining limitations.

Before production integration, extend array validation to substituent extraction,
normalization and position assignment, and floating-component localization. The
current implementation also uses recursive emission and materialized subtree
signatures (potentially quadratic string storage on long chains); this experiment
does not establish performance or robustness for adversarially deep inputs.
A public supported bridge is used for graph creation, so there is no dependency
on igraph's private C ABI or R object layout.
