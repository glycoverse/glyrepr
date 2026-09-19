# IUPAC parser batching experiment

Standalone prototype; production package code is unchanged. This experiment tests
batching in R before considering a C++ parser, and measures its interaction with
the previous Rcpp sequence kernel.

`prototype.R` retains the existing string validation, normalization, tokenization,
residue/linkage parsing, branch-stack implementation, reducing-end handling and
error wrapper. It replaces repeated per-residue igraph mutations with preallocated
residue/substituent/linkage/edge arrays, one graph allocation for all vertices, and
one edge insertion. The C++ sequence prototype is loaded from `../rcpp-prototype`.
Namespace bindings are replaced only temporarily in the running R process.

## Reproduce

Run sequentially from the repository root, without other heavy work:

```sh
Rscript benchmarks/parser-batching/run.R verify
Rscript benchmarks/parser-batching/run.R bench 1
Rscript benchmarks/parser-batching/run.R bench 2
Rscript benchmarks/parser-batching/run.R bench 3
Rscript benchmarks/parser-batching/summarize.R
```

Dependencies are the package development environment, glydb, and Rcpp plus a C++
compiler for the combined condition. `benchmarks/` is excluded from package builds.
The exact 500 ordinary IUPAC strings are reused from
`../rcpp-prototype/results/sample-iupac.txt`; this fixes input across experiments.
These timings are newly measured against their own baseline, not compared to old
elapsed times recorded under different host conditions.

Three separate scopes are timed:

- **Graph assembly:** prepared residue/linkage/endpoint arrays to a raw igraph.
  Compares incremental versus batched construction. Excludes string processing,
  branch-stack work, graph-level attributes, validation and canonicalization.
- **Complete tree parser:** IUPAC string to raw igraph, including tokenization,
  substituent/linkage parsing, branch-stack processing, assembly and reducing-end
  attributes. Compares the actual old internal parser and the batched prototype.
- **Public constructor:** `as_glycan_structure(strings)`, including parsing,
  validation, canonical ordering, IUPAC generation and vector construction. Four
  conditions: baseline, batch only, sequence kernel only, and both combined.

Each timing processes all 500 inputs. Three fresh R processes each warm every
condition then take three rounds, rotating/reversing condition order. Compilation,
preparation and explicit garbage collection happen outside the timers. Garbage
collection inside an operation is included. Speedups are within-worker/round ratios;
reported ranges are observed ranges, not confidence intervals.

Verification checks raw parser output for all 8,573 installed corpus structures,
including floating metadata. Public construction is compared across the full corpus
for batch-only and combined variants, including exact graph/vertex/edge attributes
and endpoint order. Additional checks cover names, NA, duplicates, alditols,
ambiguous linkages, substituents and floating structures. The full package suite is
run separately with batch-only and combined bindings. The isolated assembly
benchmark verifies identical graphs for every prepared input before timing.

`results/` contains raw timings, paired ratios, test results and R session details.
`provenance.json` identifies baseline source and input hashes. testthat may normalize
trailing blank lines in existing snapshots; review and remove such incidental edits.

## Results

Each elapsed time is the median for 500 inputs over nine measurements in three
fresh processes. Speedup is the median paired ratio, which can differ from the
ratio of elapsed-time medians. Ranges below retain all measurements.

| Operation / variant | Baseline median (s) | Variant median (s) | Paired speedup | Observed paired range |
|---|---:|---:|---:|---:|
| Graph assembly / batch | 0.493 | 0.113 | 4.31x | 4.21-4.93x |
| Complete parser / batch | 0.967 | 0.552 | 1.73x | 1.21-2.74x |
| Public constructor / batch only | 1.727 | 1.339 | 1.30x | 1.23-1.59x |
| Public constructor / sequence kernel only | 1.727 | 1.635 | 1.06x | 1.02-1.26x |
| Public constructor / batch + sequence kernel | 1.727 | 1.235 | 1.40x | 1.14-1.51x |

Adding the sequence kernel on top of batching gives a median paired incremental
speedup of 1.08x (observed 0.72-1.11x). The combined result is
measured directly; it is not obtained by multiplying separate speedups.

Both batching-only and combined variants passed all 1,965 package test expectations
with zero failures, errors, warnings or skips. Raw-parser parity covered 8,573 corpus
entries; public-constructor parity covered all entries for both variants. Assembly
outputs also matched for all 500 prepared benchmark inputs in each worker.

Batching in R alone is worthwhile: approximately 30% greater end-to-end throughput
(about 23% less elapsed time using the paired median), with no new compiled parser.
Combining it with the earlier Rcpp sequence kernel gives approximately 40% greater
throughput (about 29% less elapsed time). The isolated 4.3x assembly improvement
must not be presented as a 4.3x package-wide improvement. Parser timings showed
more variability than assembly timings; all raw rounds, including outliers, remain
in the results.

This supports production integration of batching before attempting a C++ scanner.
The string parser remains R code; this experiment does not measure a native C++
parser. Production package files and dependencies were not changed.
