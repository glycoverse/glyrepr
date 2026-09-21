# Public API integration and compatibility

The compact C++ parser is now a private package backend, compiled at installation.
Callers use the existing `as_glycan_structure()` and vctrs APIs; there is no new
public backend selector, argument, class or return representation.

## Routing and compatibility

| API/path | Behavior |
|---|---|
| `as_glycan_structure(character, on_failure="error")` | Parse, validate and canonicalize distinct strings using arrays; construct each distinct final graph once. |
| `as_glycan_structure(character, on_failure="na")` | The same array backend, with element-local recovery, original indices/names and one existing aggregate warning. |
| `vctrs::vec_cast(character, glycan_structure())` and character consumers | Existing S3 dispatch reaches the private array backend. |
| igraph/list conversion and `glycan_structure(...)` | Existing R graph validation/canonicalization, including arbitrary graph/vertex/edge attributes and the existing missing-input rules. |
| `new_glycan_structure()`, low-level validation/canonicalization/IUPAC helpers | Existing trusted/untrusted input contracts retained. |
| Inspection, transformations, tables, floating localization, S3/vector operations | Existing implementations consume the identical canonical graphs and keys. |

This is API compatibility coverage, **not** a claim that graph-input and every
accessor implementation were rewritten in C++. In particular, character input
was not added to `glycan_structure(...)`, and invalid list/container types retain
their existing errors. The 52 S3 registrations and all 78 public function
signatures are unchanged; NAMESPACE adds only native-library registration and an
Rcpp import. Source installation now needs Rcpp and a C++17 compiler. Installed
use does not invoke `sourceCpp()` or runtime compilation.

For strict input, if any native element fails, the complete original R path is
replayed. This preserves the original two-stage **parse all, then validate** error
precedence and the existing purrr condition indices, names and messages. Under
`on_failure="na"`, only failed native elements go through the reference parser and
validator for their original diagnostic. Existing recovery assembly determines
positions, repeats failure reasons for duplicate inputs, preserves names/NA and
emits the same warning class/message. Native size guards likewise fall back
rather than imposing new public size limits.

## Verification

Reference: commit `5b500b2`, installed in a separate temporary R library before
production edits. Candidate: the final source package built and installed by the
successful R CMD check. Comparisons ran in fresh R processes against each library.

- **78 exported signatures identical**, including defaults.
- **113 public behavior cases identical**, including result classes, values,
  graph metadata, warning/error classes, complete messages, condition calls and
  parent-condition signatures. Internal trace objects/session metadata are not
  compared.
- Covered character, igraph, graph list, existing vector, NULL, missing, empty,
  malformed, factor/numeric/matrix inputs, duplicate names/values, invalid failure
  policies and parse-before-validation precedence.
- Covered both policies, vctrs casts, graph constructors, trusted construction,
  graph access, IUPAC conversion, attributes, inspection, transformations, table
  round-trips, mapping, floating enumeration/localization and vector operations.
- **8,573 corpus entries** match the saved reference canonical keys and full
  ordered graph signatures under **both `error` and `na`**. Instrumentation of the
  installed package recorded **zero reference-parser calls** for these valid
  entries.
- Full package suite: **2,017 PASS, 0 FAIL, 0 WARN, 0 SKIP**. Existing error/warning
  snapshot content was preserved. New regression snapshots cover repeated failed
  positions and strict parsing/validation precedence.
- Final `devtools::check(document=FALSE, error_on="note")` / R CMD check
  (`--no-manual --as-cran`): **0 errors, 0 warnings, 0 notes**. Native compilation,
  loading/unloading, examples, tests and vignette rebuilding passed on arm64 macOS,
  R 4.6.1. This is a local result, not a cross-platform CI claim.

The native collation and single BLISS symmetry bridges documented in
[RESULTS-extended.md](RESULTS-extended.md) remain in place. Graph-input paths retain
their R implementation so custom attributes and existing low-level semantics are
preserved.

## Installed public-entry performance

The original frozen 500 mixed inputs were reused. The recovery workload appends
20 copies of the same invalid string and one existing NA (521 total rows).
Each method ran in three fresh processes, two repetitions per workload/process.
Reference/candidate process order alternated across worker pairs. Loading and
warmups are excluded; parsing, validation, graph creation, vector assembly and
warning construction are included. Warnings were suppressed on both sides.

| Public call/workload | R median seconds (range) | Candidate median seconds (range) | Median paired speedup |
|---|---:|---:|---:|
| 500 valid mixed, `on_failure="error"` | 1.953 (1.814–2.107) | 0.125 (0.123–0.151) | **14.98×** (13.95–16.61×) |
| 500 valid mixed, `on_failure="na"` | 1.890 (1.803–1.995) | 0.134 (0.126–0.147) | **13.88×** (13.33–14.71×) |
| 521 rows including repeated failures, `on_failure="na"` | 1.937 (1.824–2.024) | 0.178 (0.154–0.278) | **10.90×** (6.56–13.14×) |

Speedups are medians of six paired ratios, not ratios of displayed medians.
All observations, including slower runs, are retained. These are workload-specific
measurements against the separately installed old package, not against the current
public function mislabeled as an R reference.

## Evidence and reproduction

Evidence is in `results-api/`: baseline corpus signatures and public API snapshots,
[public parity](results-api/public-parity.txt), [installed corpus parity](results-api/installed-parity.txt),
[final check](results-api/check-final.log), [test output](results-api/tests-final.Rout),
[raw timings](results-api/timings.csv), [paired ratios](results-api/paired-speedups.csv),
and [source/input hashes](results-api/provenance-md5.csv). The first check's fixed
namespace NOTE and initial new-snapshot warnings are retained as historical logs;
use the final check/test artifacts for the final status.

To repeat the two-library comparison, install commit `5b500b2` and this commit to
separate temporary R libraries. Pass those library paths to the same script:

```sh
Rscript benchmarks/compact-cpp/audit-public-api.R /path/to/reference-lib /tmp/reference.rds
Rscript benchmarks/compact-cpp/audit-public-api.R /path/to/candidate-lib /tmp/candidate.rds
Rscript benchmarks/compact-cpp/audit-installed.R /path/to/candidate-lib
Rscript benchmarks/compact-cpp/benchmark-public-api.R /path/to/reference-lib 1 R
Rscript benchmarks/compact-cpp/benchmark-public-api.R /path/to/candidate-lib 1 compact
```

Repeat benchmark worker IDs 2 and 3, reversing reference/candidate process order
for worker 2, then run `summarize-public-api.R`. It checks the saved standard public
snapshots and summarizes the six paired timing files. `audit-installed.R` checks
against retained `baseline-*.rds` corpus signatures. The reference capture script
must be run on the reference checkout to regenerate those signatures; running it
on the new package would replace the oracle with the candidate.
