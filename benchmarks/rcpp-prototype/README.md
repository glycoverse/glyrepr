# Rcpp canonicalization / serialization prototype

This standalone experiment leaves the production R code, DESCRIPTION and NAMESPACE
unchanged. `benchmarks/` is already excluded from package builds. Requires the
package's development dependencies, Rcpp, a C++ compiler, and the installed glydb
corpus used in the recorded session information.

The C++ kernel accepts integer edge endpoints and residue, substituent and linkage
vectors. It builds adjacency, subtree depths and signatures, determines branch
order, and emits vertex order, edge order and IUPAC in one call. Output assembly
uses an explicit stack and a single string buffer. Subtree signatures are still
materialized, so extremely deep trees can retain quadratic signature storage.
R's `sort()` and `order()` are called at branching nodes to preserve current locale
collation and stable tie behavior; this is intentionally not a completely native
sorting implementation.

R retains graph extraction/reordering, reducing-end formatting, alditol handling,
public validation and errors, and floating metadata. Graphs with floating metadata
attributes use the original implementation. The prototype's C++ entry point is
internal and assumes validated glycan labels/linkages; it is not a public parser.

`with_prototype()` temporarily replaces two namespace bindings within one process
and restores them on exit. This is a benchmark harness, not a production integration
mechanism. Compilation uses a temporary directory.

## Reproduce

From the package root:

```sh
Rscript benchmarks/rcpp-prototype/run.R verify
Rscript benchmarks/rcpp-prototype/run.R bench 1
Rscript benchmarks/rcpp-prototype/run.R bench 2
Rscript benchmarks/rcpp-prototype/run.R bench 3
Rscript benchmarks/rcpp-prototype/summarize.R
```

Run these sequentially without other CPU-heavy work. Validation compares sequence
orders, strings, graph/vertex/edge attributes and edge endpoints across the full
installed corpus; it also checks 500 seeded vertex/edge permutations, targeted
fixtures, three available collation locales, names/NA/duplicates, and the complete
package test suite with the prototype bindings enabled. testthat may normalize
trailing blank lines in existing snapshot files; review any working-tree changes.

The timing input is 500 ordinary graphs sampled with seed 123, plus seeded
permutations and their canonical IUPAC strings. `results/sample-iupac.txt` records
the exact strings. Three fresh R processes each warm both implementations and run
three rounds, alternating implementation order. Compilation, corpus preparation,
and explicit pre-run garbage collection are outside the timers; garbage collection
during an operation is included. Each operation processes all 500 graphs/strings.
The sequence-kernel measurement includes R-side graph extraction for both paths.

`results/timings-*.csv` contains all raw elapsed times; `paired-timings.csv` and
`summary.csv` contain paired speedups and median elapsed times. These are descriptive
local benchmarks, not confidence intervals. `session-*.txt` records software and
locale, `validation.txt` records parity/test counts, and `tests.csv` records each
package test result. Floating workloads receive fallback coverage but have no
claimed acceleration. Source-string parsing and graph reconstruction remain in R,
so their end-to-end improvement is expected to be smaller than serialization's.

## Recorded results

All values are for 500 structures. There are nine paired measurements per operation
across three fresh processes. Speedup is the median of within-round R/Rcpp ratios,
not the ratio of the two elapsed-time medians; these differ under timing drift.

| Operation | R median (s) | Rcpp median (s) | Paired speedup median | Observed paired range |
|---|---:|---:|---:|---:|
| Construct from canonical graphs | 0.479 | 0.381 | 1.26x | 1.20-1.37x |
| Graph to IUPAC | 0.179 | 0.087 | 2.08x | 1.94-2.17x |
| Construct from IUPAC strings | 2.231 | 2.253 | 1.05x | 0.92-1.07x |
| Construct from permuted graphs | 0.956 | 0.798 | 1.10x | 0.90-1.24x |
| Sequence kernel, including extraction | 0.171 | 0.076 | 2.23x | 2.18-3.00x |

Validation: 8,009 ordinary graphs matched kernel vertex/edge orders and strings;
all 8,573 graphs matched canonical graph attributes, endpoints and IUPAC, including
564 floating-metadata fallbacks. Another 500 vertex/edge permutations passed.
Targeted fixtures and 50 permuted graphs passed under C.UTF-8, C and en_US.UTF-8.
The prototype-enabled package suite passed 1,965 expectations, with zero failures,
errors, warnings or skips. Production source files remain unchanged.

Interpretation: serialization and canonical-graph construction show consistent
improvements in these samples. Permuted-graph construction is variable, including
some slower pairs. String construction has no convincing end-to-end improvement:
its paired median is only 1.05x, some pairs regress, and the aggregate elapsed-time
median is slightly worse. Host timing drift was visible, particularly in the last
worker; the raw pairs are retained rather than discarding unfavorable observations.

Recommendation: this prototype supports integrating an ordinary-tree sequence
kernel for serialization and canonicalization, with the existing floating fallback.
Do not claim a general 2x package speedup. Parser batching and graph rebuilding are
separate opportunities. Before production integration, replace the temporary
namespace-binding harness with registered package Rcpp bindings and dedicated
regression tests, and run package/platform and reverse-dependency checks.
