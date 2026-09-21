# Results: compact C++ IUPAC pipeline

Measured 2026-09-21 against checkout `0b7231bac1b5ed207ee9e55524f928534233ab76`
(the existing batched pure-R parser). The new design delivers a substantial
end-to-end improvement on supported trees, with smaller gains once unsupported
syntax invokes the existing pipeline. No production functions or dependencies
were changed.

## End-to-end constructor

Each row has six paired measurements across three fresh, sequential R processes.
Times are median elapsed seconds with observed min–max. Speedup is the **median
of paired R/compact ratios**, so it need not equal the ratio of the two medians.

| Workload | Current R, seconds | Compact, seconds | Paired speedup (range) |
|---|---:|---:|---:|
| 500 supported ordinary trees | 1.235 (1.198–2.494) | 0.121 (0.110–0.240) | **10.77×** (9.81–14.48×) |
| 500 ordinary trees, including 44 modified-residue fallbacks | 1.699 (1.258–2.179) | 0.329 (0.260–0.695) | **4.92×** (3.14–5.32×) |
| 500 mixed entries: 433 fast, 32 floating, 35 modified | 2.204 (1.971–3.766) | 1.144 (0.925–1.707) | **2.01×** (1.89–3.05×) |
| 10,000 rows repeating 100 supported inputs | 0.292 (0.274–0.529) | 0.025 (0.024–0.044) | **11.90×** (11.40–15.77×) |

These timings include the complete strict character-input constructor, including
validation, graph creation, canonical keys, deduplication, vector assembly and
fallbacks. The duplicate workload benefits from within-call deduplication in
**both** implementations; no parsed-input cache persists between measurements.
Ordinary/mixed samples were drawn independently with seed 20260921 and retained
as snapshots. These are sample-specific results, not full-corpus timing claims.

## Isolated stages on 500 supported trees

| Stage | Current R median (range), seconds | Compact median (range), seconds |
|---|---:|---:|
| Parse + validate + canonical graph/IUPAC pairs | 2.043 (1.460–3.759) | 0.145 (0.119–0.212) |
| Array parsing + validation + canonical order/IUPAC only | — | 0.022 (0.019–0.033) |
| Materialize precomputed arrays as final igraphs | — | 0.147 (0.104–0.218) |

The equivalent graph/IUPAC-pair stage is **13.12×** faster by paired median
(range 9.31–29.37×). The separate array and materialization stages produce
different outputs; their timings are absolute costs, not a speedup comparison.
They are measured separately and should not be added to reconstruct another
stage's median. GC and host variability are evident, especially in worker 1 and
the slow R pair-stage measurement; all measurements are retained.

A separate 1 ms sampled profile attributed 80.1% of compact-constructor time to
`compact_graph()`, including 61.0% in `igraph::make_graph()`, and 15.8% to the
array kernel. Percentages are inclusive sampled time from one profiled run,
not confidence intervals or an allocation measurement. Graph materialization
is now the main remaining supported-path cost.

The instrumented audit on 1,001 input rows (500 supported inputs repeated twice,
plus NA) observed exactly **500 `make_graph()` calls**, and **zero** calls to
`make_empty_graph()`, `add_edges()`, `permute()`, or `graph_from_data_frame()`.
Thus each distinct supported canonical graph is constructed once, already in
final order; attributes are assigned in batches afterward.

## Correctness and scope

- All **8,573** glydb entries matched the reference canonical string and exact
  ordered graph endpoints, vertex attributes, edge attributes and graph metadata.
- **7,356 (85.8%)** used C++ arrays; **564** floating entries and **653** entries
  with modified residues used the reference fallback. Full-corpus parity for
  fallback entries does not imply C++ implementation of those features.
- Every known residue and its alditol, missing reducing-end annotations,
  ambiguous/unknown linkages, names, NA, empty vectors, raw duplicates and
  canonical-equivalent branch permutations were checked.
- **500** seeded random trees of 2–40 residues matched, including noncanonical
  branches and signature ties. They also matched under `C`, `C.UTF-8` and
  `en_US.UTF-8` collation.
- **17** malformed inputs matched reference acceptance/rejection. Of **304**
  string literals extracted from existing parser/serializer tests, all **61**
  fast-accepted literals were accepted by the reference with identical keys.
- No package-wide test/check claim is made: package code is unchanged, and the
  experiment is verified by its standalone parity and instrumentation scripts.

The fast path validates its bounded grammar on arrays; it does not merely trust
corpus strings. The non-C collation bridge still invokes `base::order()` for
signature ordering. Modified residues, floating metadata, graph input,
`on_failure = "na"`, and deep/adversarial-input engineering remain outside the
native implementation. Those limits should be addressed before integration.

## Environment and evidence

R 4.6.1, Rcpp 1.1.2, igraph 2.3.3, arm64 macOS Tahoe 26.5.2; timing locale
`C.UTF-8`. Fresh-process source/compilation took 5.285–7.252 seconds, recorded
separately and excluded from steady-state measurements. First-use runtime
compilation could therefore dominate a small one-off call.

Reproduction commands and design are in [README.md](README.md). Evidence includes
[raw timings](results/timings.csv), [paired ratios](results/paired-speedups.csv),
[coverage](results/coverage.csv), [parity](results/parity.txt),
[graph-call counts](results/graph-calls.csv), and
[source/input checksums](results/provenance-md5.csv). Corpus, random inputs,
workloads, profiles and per-process session details are retained alongside them.
