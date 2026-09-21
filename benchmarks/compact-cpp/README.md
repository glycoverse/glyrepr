# Current package integration

The array backend is now integrated behind the unchanged public APIs, including
`on_failure = "na"`. See [API-COMPATIBILITY.md](API-COMPATIBILITY.md) for compatibility,
installed-package verification, dependencies and public-entry benchmarks.

The material below describes the historical standalone prototype. Reproduce its
v2 experiment at commit `5b500b2` (v1 at `e735332`); on the current checkout,
`as_glycan_structure()` is already accelerated and is **not** a pure-R benchmark
reference. Current comparisons use separately installed reference/candidate packages.

---

# Compact-array C++ IUPAC prototype

The current experiment parses ordinary trees, modified residues, floating glycan
parts and floating substituents into C++ arrays. It validates and canonicalizes
those arrays before creating each distinct final glycan igraph once. Package
code and dependencies remain unchanged.

The first experiment is preserved in [RESULTS.md](RESULTS.md) and `results/`.
It covered 7,356/8,573 inputs natively and fell back for the other 1,217.
The extension is described in [RESULTS-extended.md](RESULTS-extended.md), with
fresh evidence in `results-v2/`. The old version is reproducible at commit
`e735332`; the original benchmark scripts describe that historical experiment.

## Current implementation

`compact.cpp` owns residue/substituent/linkage/parent arrays, adjacency arrays,
subtree signatures, canonical traversal and floating metadata. It handles:

- Known concrete/generic residues, unusual configurations, reducing-end inference
  and alditols; ambiguous and unknown linkage positions.
- Substituent extraction, Neu5Ac/Neu5Gc rules, normalized slash-separated
  positions, stable ordering and conflict-free substituent position assignment.
- Leading floating blocks and global input parent indices, including floating
  substituents, parent domains across components, and explicit parent validation.
- Joint carbon-slot feasibility via bipartite matching, component reachability,
  acyclic parent assignment, and iterative singleton-component attachment.
- Canonical component/branch ordering, symmetry-aware ties, complete-sequence
  parent remapping, and final graph metadata and IUPAC serialization.

No glycan igraph is created during parsing or validation. The final R bridge uses
one `igraph::make_graph()` per distinct canonical glycan and attaches attributes
in batches. It never inserts edges, permutes a graph, or reconstructs one.
Duplicates and NA preserve vector semantics, and there is no cross-call input
cache. Vocabulary/configuration tables come from the package, not the corpus.

Two explicit library boundaries remain. Non-C collation uses `base::order()` on
small signature vectors, as in the original prototype. A forest with unresolved
floating metadata makes **one BLISS callback** through public igraph APIs. C++
builds the colored augmented-constraint arrays; the callback constructs one
auxiliary undirected graph and obtains canonical labels. C++ then completes
ordering and metadata remapping. This is not a claim of zero R callbacks or one
graph allocation including the auxiliary graph. It avoids per-component glycan
graph creation and all repeated R parsing/validation/canonicalization work, and
does not depend on igraph's private object layout or C ABI.

The character constructor's `fallback = FALSE` mode strictly exercises this
array pipeline. Verification and benchmarks use that mode. Default fallback
exists only to get the public error diagnostic for a native failure or size
limit; passing benchmarks must never use it to hide unsupported valid inputs.
The prototype has a 2,048-residue limit per complete structure. Recursive
traversal, materialized subtree signatures and potentially combinatorial
floating-parent search still need deep/adversarial-input engineering before
production. It does not implement graph-input conversion or `on_failure="na"`.

## Reproduce the extension

Run from the package root with package development dependencies, Rcpp and glydb.
Use sequential workers to avoid benchmark contention:

```sh
Rscript -e 'source("benchmarks/compact-cpp/verify-extended.R")'
Rscript -e 'source("benchmarks/compact-cpp/audit-extended.R")'
Rscript -e 'source("benchmarks/compact-cpp/audit-graph-calls-extended.R")'
Rscript benchmarks/compact-cpp/benchmark-extended.R 1
Rscript benchmarks/compact-cpp/benchmark-extended.R 2
Rscript benchmarks/compact-cpp/benchmark-extended.R 3
Rscript benchmarks/compact-cpp/audit-size-guard.R
Rscript benchmarks/compact-cpp/summarize-extended.R
```

Each benchmark worker starts in a fresh R process and records source/compilation
time separately. Six paired repetitions alternate baseline/prototype order;
GC and three-element warmups occur outside timing regions. Four original
workloads are reused verbatim. New workloads add 500 modified inputs selected
from the original 653 fallback entries with seed 20260922, and all 564 floating
entries. All timings include the native collation/BLISS bridges as applicable.

The constructor comparisons include parsing, validation, ordering, graph creation,
canonical deduplication and vector assembly. The pair stage returns equivalent
canonical graph/IUPAC pairs on both paths without vector assembly. Array-only
and final-materialization timings are separate absolute measurements, not
comparisons of equivalent outputs or additive decompositions of medians.

Evidence is retained as corpus/workload/random-input RDS snapshots, per-input
coverage/parity CSVs, graph-call counts, six raw timing repetitions, paired
ratios, session details, profiles and source/input checksums. The extended audit
mines existing test literals and adds modified-residue, floating-parent,
cycle/conflict, singleton and symmetry cases; it compares native acceptance and
rejection directly to R without fallback.
