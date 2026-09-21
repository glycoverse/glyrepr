# Array backend validation and performance

Source implementation: `9b4db7ee64e4375ea7d5882bc5fa0a85c6a310d4`.
Raw timings, source hashes, session information, workload, and per-input parity
are in `results/`. Timings use an installed release build on this machine,
not the development/debug DLL. They are workload-specific, not a claim about
whole-format parsing speed in glyparse.

## Correctness

- 11,425 distinct frozen strings inspected.
- 10,640 R-reference-valid records passed after reversing source node IDs and
  remapping floating metadata. All 10,640 had native status `ok`; no valid record
  needed graph recovery. Canonical keys, vertex and edge attributes, endpoints,
  reducing-end state, and floating metadata matched the reference.
- 785 reference-invalid strings were recorded separately.
- Final glyrepr test suite: 2,084 passes, no failures or warnings.
- glyparse tests against the installed glyrepr build: 1,665 passes, no failures
  or warnings. Tests ran in a temporary copy; glyparse source was not changed.
- Installed-session public-export and named/missing round-trip smoke passed.
- Documentation/reference-index checks passed.
- Final `R CMD check --no-manual`: 0 errors, 0 warnings, 2 environment notes
  (remote clock verification and the toolchain's temporary `xcrun_db`). The
  devtools wrapper also reported a failed `quarto -V` probe; package examples,
  tests, and vignette checks completed successfully.

The two pre-existing canonical-graph revalidation failures are documented in
[README.md](README.md#existing-reference-round-trip-limitation), with retained
reproductions under `known-reference/`. They are not hidden by the successful
raw-record parity result.

## Complete construction from parsed records

All times are seconds. Each cell is median [minimum, maximum] across three
alternating-order rounds after warmup. Format-specific string parsing is outside
both measured constructors; both produce complete `glyrepr_structure` vectors.

| Workload | Build graphs + existing constructor | New array constructor | Ratio of medians |
|---|---:|---:|---:|
| 500 mixed records, 31 with floating metadata | 11.325 [10.943, 15.084] | 0.257 [0.210, 0.277] | 44.1x |
| 1,000 records repeating 20 source records | 1.851 [1.760, 1.861] | 0.151 [0.149, 0.153] | 12.3x |

Paired per-round speedups for the mixed workload range from 44.1x to 54.5x.
The repeated workload still validates and canonicalizes every array record;
canonical graph materialization/storage is deduplicated. glyparse can also
continue deduplicating original strings before parsing.

## Existing graph callers

| Workload | Separate validation/order/key calls | Fused graph batch |
|---|---:|---:|
| 500 mixed records | 14.966 [10.772, 15.077] | 14.959 [10.787, 16.324] |
| 1,000 repeated records | 1.502 [1.430, 1.531] | 1.526 [1.481, 1.585] |

These measurements do not establish a material end-to-end speedup for the graph
batch API. It consolidates the interface and avoids duplicate ordinary-tree
traversal while retaining arbitrary attributes; floating normalization still
uses the reference graph path. Direct array output is the main performance
migration for glyparse.

Native-array-only medians are 0.019 and 0.011 seconds respectively. They exclude
R record validation and igraph materialization and are not complete-constructor
speedups.
