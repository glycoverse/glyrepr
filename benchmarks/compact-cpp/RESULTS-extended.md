# Expanded compact-array prototype: modified and floating structures

Measured 2026-09-21 against the unchanged package R implementation at checkout
`e735332522feeaca2f4009f78dd6b8e113400e6a`. The extended array pipeline now covers
**all 8,573 corpus entries**, including the 653 modified-residue entries and
564 floating entries that previously used R parsing fallback. Package production
code, DESCRIPTION and NAMESPACE remain unchanged.

## Correctness and coverage

Full-corpus verification ran with `fallback = FALSE`. Every canonical string,
ordered edge endpoint, vertex attribute (including `sub`), edge attribute,
graph attribute and floating metadata field matched the reference. Coverage is
8,573 `ok`, zero fallback and zero errors. Names, NA, empty vectors, duplicates,
canonical-equivalent branch permutations, all known residues and their alditols,
and 500 seeded noncanonical random trees also passed. Random-tree comparisons
ran under C, C.UTF-8 and en_US.UTF-8 collation.

The independent extended audit tested **2,870 unique inputs**, combining existing
parser/floating/substituent/alditol test literals, residue/substituent combinations,
400 generated floating cases, and targeted constraints. Native acceptance matched
R for **2,085 accepted** and **785 rejected** inputs. Every accepted result had
identical canonical strings and full ordered graph signatures. Valid targeted
floating symmetry cases also passed under the three collation settings. This
audit calls `compact_arrays()` directly and cannot conceal native errors behind
R fallback. The existing 17 malformed cases were also retained in the corpus
verification script. A separate boundary check confirms the complete-forest
2,048-residue guard.

Covered constraint cases include ambiguous substituent positions, repeated unknown
positions, Neu5Ac/Neu5Gc parsing, carbon-position conflicts, explicit impossible
parents, component-to-component attachments, closed cycles, chained singleton
resolution, floating substituents, symmetry ties, and global parent-index remapping.

## End-to-end performance

Six paired repetitions across three fresh sequential R processes. Values are
median elapsed seconds with min–max. Speedups are **medians of paired ratios**
with their observed min–max, not ratios of the displayed marginal medians.
Compilation is excluded. Every prototype measurement disables R parsing fallback.

| Workload | Current R seconds | Extended prototype seconds | Paired speedup |
|---|---:|---:|---:|
| 500 modified-residue inputs | 1.865 (1.824–1.999) | 0.112 (0.111–0.158) | **16.61×** (11.54–18.01×) |
| All 564 floating inputs | 22.022 (21.979–26.115) | 0.252 (0.247–0.389) | **87.79×** (67.13–89.10×) |
| Original 500 mixed inputs | 1.981 (1.837–2.071) | 0.128 (0.123–0.206) | **15.55×** (9.30–16.43×) |
| Original 500 ordinary inputs | 1.332 (1.306–1.370) | 0.116 (0.115–0.117) | **11.49×** (11.35–11.71×) |
| Original 500 unmodified supported inputs | 1.271 (1.239–1.364) | 0.116 (0.115–0.124) | **11.00×** (9.99–11.75×) |
| Original 10,000 rows / 100 distinct inputs | 0.249 (0.247–0.258) | 0.023 (0.022–0.023) | **10.85×** (10.74–11.23×) |

The original four workloads are reused from their saved RDS snapshot, avoiding
selection drift. The modified group samples 500 of the original 653 entries with
seed 20260922. The floating group includes all 564 original floating entries.
The original mixed workload has 433 originally supported, 35 modified, and 32
floating inputs; all now use the array path. These are comparisons against a
freshly timed R baseline, not direct cross-session timing ratios against v1.

Timings include parsing, validation, canonical ordering, symmetry callbacks,
final graph creation, deduplication and structure-vector assembly. Both paths
deduplicate repeated input strings within a call. There is no persistent parsed
input cache. Retained high measurements show host/GC variability; no run was
removed and no package-wide performance claim is implied.

## Isolated stages and graph calls

On the original 500 unmodified inputs, equivalent parse/validate/canonical
**graph-IUPAC pairs** took 1.304 s (1.288–1.893) in R versus 0.116 s
(0.112–0.157) in the prototype: paired median **11.62×** (10.78–12.06×).
Array parsing/validation/order/IUPAC alone took 0.024 s (0.023–0.032), and final
graph materialization from those arrays took 0.090 s (0.088–0.128). These last
two stages have different outputs and are absolute timings, not a speedup
comparison or an additive decomposition of medians.

An instrumented run used all 1,217 formerly unsupported inputs plus 50 ordinary
inputs, repeated twice, plus NA: 2,535 rows, 1,267 distinct final graphs.

| Operation | Calls |
|---|---:|
| Final glycan graph creation | **1,267** |
| Auxiliary BLISS graph creation | **564** |
| `igraph::make_graph()` total | 1,831 |
| `igraph::canonical_permutation()` | 564 |
| `make_empty_graph()`, `add_edges()`, `permute()`, `graph_from_data_frame()` | **0 each** |

Thus the final glycan is still constructed only once, already in final order.
There are no temporary per-component glycan graphs, graph permutations, or graph
reconstructions. An unresolved floating forest additionally creates one auxiliary
colored constraint graph for canonical symmetry labels.

## Implementation boundaries

Residue/substituent extraction, carbon-domain matching, floating-parent
reachability/acyclicity, singleton localization, branch/component ordering,
parent remapping and output serialization operate on C++ arrays. Validation is
performed on the supplied strings and constraints, not trusted corpus keys.

**This is not an all-native, zero-callback backend.** Non-C string collation still
uses `base::order()` on signature vectors. Unresolved floating metadata uses one
public igraph BLISS callback; its auxiliary graph creation and labeling costs are
included in all relevant timings. The callback does not call glyrepr parsing,
validation or canonicalization routines. C++ supplies the augmented graph arrays
and consumes the labels to finish canonical ordering. Neither bridge accesses
igraph's private representation or ABI.

The default wrapper can still invoke the reference solely for a native failure's
public diagnostic or a size-limit case. The reported accepted-input audits and
benchmarks disable this behavior. The prototype supports strict character input,
not graph-input conversion or `on_failure = "na"`. Deep trees, potentially
quadratic subtree-signature storage and combinatorial parent-domain search still
need production engineering; matching itself uses augmenting paths rather than
exhaustive carbon-slot enumeration. Existing package validation semantics are
preserved rather than adding new biochemical restrictions.

## Reproduction and evidence

R 4.6.1, Rcpp 1.1.2, igraph 2.3.3, arm64 macOS Tahoe 26.5.2, C.UTF-8 timing
collation. Fresh source/compilation took 4.312–7.513 s and is excluded from the
steady-state results. Initialization can therefore dominate small one-off calls.

See [README.md](README.md) for commands. Retained evidence:

- [Corpus parity](results-v2/parity.txt) and [per-row coverage](results-v2/coverage.csv).
- [Extended audit](results-v2/extended-audit.txt), [per-input results](results-v2/extended-audit.csv), and frozen input RDS files.
- [Raw timings](results-v2/timings.csv), [paired ratios](results-v2/paired-speedups.csv), and [summary](results-v2/summary.csv).
- [Graph-call audit](results-v2/graph-calls.csv) and [size guard](results-v2/size-guard.txt).
- Per-worker sessions, sampled profiles, saved workloads, and [source/input hashes](results-v2/provenance-md5.csv).

Only standalone prototype checks are claimed; no package-wide test/check was
needed or run because production package code is unchanged.
