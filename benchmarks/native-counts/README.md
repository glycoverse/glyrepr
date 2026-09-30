# Native composition and count audit

Run from the repository root with the newly installed package on `.libPaths()`:

```sh
R_LIBS=/path/to/installed/library Rscript benchmarks/native-counts/audit.R
```

The comparison uses the retained R graph-to-composition implementation plus
`smap()` and the unchanged composition method of `count_mono()`. Both paths run
in the same installed package with identical dependencies. This isolates this
change from earlier package improvements.

Inputs cover 300 entries containing three repeated structures and 50 distinct
linear structures of increasing length. Each operation checks exact equality
before timing. Five samples use adaptively chosen repetitions (up to 100) to
reduce timer quantization. Results are workload-specific; repeated inputs benefit
from counting each unique graph once instead of constructing and processing a
composition at every vector position.

`results.csv` contains median timings and speed ratios. `audit.rds` retains the
inputs, individual timing samples, repetition counts, and session information.
The public API and composition-input methods are unchanged. Unusual low-level
graph representations retain the existing R fallback.
