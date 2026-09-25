# Structure progress validation

Baseline: commit `41310c8`, installed from a clean archive with `R CMD INSTALL --preclean`.
Current and baseline run in separate R sessions on the same macOS machine.
Corpus: the 200 structures in `benchmarks/native-graph/results/corpus.rds`.
Character and array cases contain 804 positions (four repeats plus missing values);
graph cases contain 200 structures. Each timing is the median of five runs.
CLI rendering is suppressed with `cli.progress_show_after = Inf`, while progress
bookkeeping remains enabled. These measurements do not include terminal rendering.

| Input | Failure policy | Baseline (s) | Current, off (s) | Current, on (s) |
|---|---|---:|---:|---:|
| strings | error | 0.157 | 0.132 | 0.124 |
| strings | na | 0.127 | 0.140 | 0.110 |
| arrays | error | 1.056 | 1.052 | 1.039 |
| arrays | na | 0.839 | 1.064 | 1.018 |
| graphs | error | 0.342 | 0.309 | 0.294 |
| graphs | na | 0.380 | 0.388 | 0.376 |

Timings remain noisy: enabled runs sometimes look faster despite additional work.
Do not interpret this as a speedup or as a tight upper bound on overhead. Disabled
runs range from 0.84 to 1.27 times baseline in this small workload. Batch boundaries,
unique-input parsing, and canonical graph caching are unchanged; disabled native
reporters make no R callbacks or clock queries.

All 12 output signatures (values, names, vertex/edge/graph attributes, topology)
match the frozen baseline exactly. Both failure policies are covered.

Validation:

- Full test suite: 2264 passed, no failures, warnings, or skips.
- Final progress-focused suite: 46 passed, including recovery diagnostics,
  fallback, callback failures, and cleanup on success, error, and interrupt.
- `pkgdown::check_pkgdown()` passed.
- `R CMD check --no-manual`: zero errors and warnings; environment notes for
  unavailable remote clock verification and temporary `xcrun_db`.
- Existing error-call attribution was also checked after the final small fix.

Reproduce from the repository root:

```sh
Rscript benchmarks/progress/audit.R /path/to/baseline/library baseline /tmp/progress-audit
Rscript benchmarks/progress/audit.R /path/to/current/library current /tmp/progress-audit
```

The current run automatically compares all output signatures with the saved baseline.
