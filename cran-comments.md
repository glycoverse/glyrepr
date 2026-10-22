## Test environments

* GitHub Actions: Linux R-devel, R-release and oldrel-1 through oldrel-4; Windows and macOS R-release.
* Cloud Linux, R 4.6.1: full remote and PDF/HTML manual checks.
* win-builder, R-devel 2026-10-21 r90733 ucrt: full checks including PDF/HTML manual.

## R CMD check results

0 errors | 0 warnings | 1 note

The remote incoming check notes 7 updates in the past 6 months. This planned minor release adds optional progress and recovery-warning controls, improves parsing and structure-operation performance, and fixes a pre-existing intermittent GC protection problem when parsing floating structures.

## Reverse dependency checks

We checked all 3 current CRAN reverse dependencies against glyrepr 1.1.0 and 1.2.0: glydraw 0.9.0, glymotif 1.0.1, and glyparse 0.8.1.

All six R CMD checks completed with 0 errors, 0 warnings, and 0 notes. We saw no new problems and failed to check no packages.
