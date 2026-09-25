# Native graph pipeline audit

The baseline is commit `59fbdc3` (before this implementation). The native implementation starts at `d12112e`; the final audited implementation is `c9ee177`. Both packages use version `1.1.0.9000`; use separate library directories to avoid accidentally comparing the same installation.

## Reproduce

Install the baseline checkout into one library and this checkout into another, then run from this checkout:

```sh
Rscript benchmarks/native-graph/audit.R baseline /path/to/baseline-library
Rscript benchmarks/native-graph/audit.R current /path/to/current-library
```

The audit uses 200 distinct valid structures extracted from the package's test literals. It checks 18 workloads, exported function signatures, representative strict/recovery conditions, and all corpus floating examples with at most 64 raw candidate combinations. Examples exceeding that limit compare the resulting error instead of enumerating them. Graph comparisons include vertices, edges, custom source attributes, graph metadata, names, missing values, and deduplication. The graph-only localization comparisons retain source vertex IDs and assignment order.

Each timing is the median of five measurements in a fresh R process with the indicated installed library. Short floating workloads run ten times per measurement and report time per call. Run baseline and current sequentially, without concurrent checks or builds. These are workload-specific wall-clock measurements, not a universal speedup guarantee.

## Implementation boundaries

- Graph access uses public igraph APIs; no internal igraph layout is read.
- Validation, canonicalization, and serialization have separate native modes. Validation returns the original graph; serialization does not resolve singleton metadata or reorder parent domains.
- Native vertex/edge mappings restore arbitrary source attributes. Newly localized edges receive the same missing attributes as the reference graph operations.
- Four built-in structure transforms batch only changed unique graphs. The structure mapper family shares native result rebuilding, while user callbacks still run in R.
- Table rows become arrays before graph construction. Native enumeration rejects incompatible assignments before materializing graphs; partial localization preserves existing vertex and edge IDs.
- The R reference paths remain for established diagnostics and unsupported representations, including graphs over the native size guard, attribute-bearing standard vectors, and unusual linkage representations. Existing array/parser behavior is retained.
- Floating symmetry continues to use the public igraph BLISS callback, and non-C collation continues to use R's ordering callback.

## Evidence

`results/` contains the corpus, baseline/current output records, five-measurement summaries, session details, and the parity diff (empty when every comparison passes). The numerical comparison is recorded in `results/comparison.csv` after the final run.

## Validation

The final implementation passes 2,218 test expectations and `R CMD check` with zero errors, warnings, or notes. The check disables remote system-clock verification and excludes only Apple's `xcrun_db` compiler cache from the temporary-directory detritus scan. All code, compilation, examples, tests, documentation, and vignette checks remain enabled. The default sandboxed run also passed those checks but reported the remote-clock and compiler-cache environment notes.

The independent installed-package audit also compares all 80 exported interfaces and 51 floating examples. During development it detected and led to a regression fix for input names on recovery-warning positions. See the saved check log, test output, and `manifest.json` for the exact checked code revision and environment settings.

## Final installed-package measurements

| Workload | Baseline (ms) | Native (ms) | Speedup |
|---|---:|---:|---:|
| `generic` | 380.0 | 98.0 | 3.88x |
| `remove_linkages` | 544.0 | 106.0 | 5.13x |
| `remove_substituents` | 106.0 | 27.0 | 3.93x |
| `fill_anomer_pos` | 622.0 | 87.0 | 7.15x |
| `graph_constructor` | 899.0 | 219.0 | 4.11x |
| `batch_canonicalize` | 1040.0 | 215.0 | 4.84x |
| `canonicalize` | 697.0 | 98.0 | 7.11x |
| `validate` | 435.0 | 133.0 | 3.27x |
| `serialize` | 121.0 | 27.0 | 4.48x |
| `structure_to_iupac_graph` | 1287.0 | 149.0 | 8.64x |
| `tables` | 1273.0 | 211.0 | 6.03x |
| `smap` | 876.0 | 155.0 | 5.65x |
| `smap2` | 1103.0 | 176.0 | 6.27x |
| `spmap` | 1150.0 | 174.0 | 6.61x |
| `simap` | 1121.0 | 178.0 | 6.30x |
| `localize` | 13.8 | 2.8 | 4.93x |
| `enumerate` | 132.3 | 17.2 | 7.69x |
| `enumerate_graph` | 101.6 | 13.0 | 7.82x |
