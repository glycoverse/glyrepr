# Structure construction benchmark

Measured on 2026-09-08 with R 4.6.1, igraph 2.3.3 and glydb 0.6.0.9000.
The baseline is `f2323233a2e44d659a55389ee2000393c571936b`; the optimized
implementation is `eef7cab`.

The input contains 500 ordinary structures sampled with seed 20260908 from
8,009 distinct ordinary `glydb::glydb_structures()` entries. The same serialized
input is reused for both versions. Raw graphs retain parser vertex order;
canonical graphs have already passed through the baseline constructor. Repeated
inputs contain 10 structures repeated 100 times.

Each version ran in a fresh R process for each of three paired comparisons,
with order baseline/optimized, optimized/baseline, baseline/optimized. Timing
excludes loading the package, preparing input, garbage collection, and checking
output equivalence. The table gives median elapsed seconds and the median of
paired speedups. Raw measurements are in `structure-construction-results.csv`.

| Input | Baseline (s) | Optimized (s) | Paired speedup |
|---|---:|---:|---:|
| 500 distinct IUPAC strings | 3.039 | 1.651 | 1.84x |
| 500 raw graphs | 0.911 | 0.694 | 1.31x |
| 500 canonical graphs | 0.624 | 0.443 | 1.42x |
| 1,000 repeated strings | 0.053 | 0.030 | 1.80x |
| 1,000 repeated graphs | 1.842 | 1.369 | 1.34x |
| 500 strings with `on_failure = "na"` | 3.234 | 1.767 | 1.86x |

Constructing 10,000 missing elements separately isolates container overhead.
Median cumulative R allocation measured with `Rprofmem()` fell from 573.69 MiB
to 0.54 MiB. This is allocation volume, not peak resident memory or the cost of
constructing valid glycans. The trusted constructor was also measured, but its
elapsed times were below the timer resolution, so no speedup is reported.

All baseline/optimized output signatures matched in every pair, including
IUPAC strings, names, vertex and edge tables, graph attributes, and deduplicated
graph storage. Existing package tests separately cover floating structures,
substituents, alditols, malformed input and warning/error recovery. These timings
do not establish performance for complex floating structures.

Validation: 1,967 test expectations passed with no failures, warnings or skips;
`R CMD check` completed with zero errors, warnings and notes. No graph validation
was removed, and the graph storage representation is unchanged.

## Reproduce

Run from the repository root. The script needs `devtools`; input preparation
also needs `glydb`. Use one input file for all runs, and run benchmarks
sequentially without competing R jobs.

```sh
mkdir -p /tmp/glyrepr-before /tmp/glyrepr-bench
git archive f2323233a2e44d659a55389ee2000393c571936b | tar -x -C /tmp/glyrepr-before
Rscript benchmarks/structure-construction.R prepare /tmp/glyrepr-before /tmp/glyrepr-bench/input.rds
Rscript benchmarks/structure-construction.R measure /tmp/glyrepr-before /tmp/glyrepr-bench/input.rds /tmp/glyrepr-bench/baseline.rds
Rscript benchmarks/structure-construction.R measure . /tmp/glyrepr-bench/input.rds /tmp/glyrepr-bench/optimized.rds
```

Repeat the last two commands in alternating order, retaining each output under a
separate filename. Compare a pair in R:

```r
baseline <- readRDS("/tmp/glyrepr-bench/baseline.rds")
optimized <- readRDS("/tmp/glyrepr-bench/optimized.rds")
stopifnot(identical(baseline$signatures, optimized$signatures))
merge(baseline$timings, optimized$timings, by = "case",
      suffixes = c("_baseline", "_optimized"))
c(baseline = baseline$missing_allocated_mib,
  optimized = optimized$missing_allocated_mib)
```

Each result also records `sessionInfo()`. Input membership may change with the
installed glydb version; use the saved input file when comparing implementations.
