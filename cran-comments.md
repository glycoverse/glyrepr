## Test environments

* GitHub Actions: Linux R-devel, R-release and oldrel-1 through oldrel-4; Windows and macOS R-release.
* macOS arm64, R 4.6.1: full remote and PDF/HTML manual checks.
* win-builder, R-devel 2026-11-07 r90850 ucrt: full checks including tests, vignettes and PDF/HTML manual, completed 2026-11-10.

## R CMD check results

0 errors | 0 warnings | 1 note

Pre-release checks used version 1.2.0.9000. The incoming NOTE identifies the development version number and 7 updates in the past 6 months. The submission version is 1.2.1.

This patch fixes a reproducible Windows package-loading hang with Chinese LC_COLLATE settings. It replaces locale-dependent native regex parsing with ASCII grammar matching without changing the session locale. It also fixes incorrect result reuse for colliding structure-operation arguments and declares the minimum supported igraph version. The Windows loading failure motivates this patch despite the recent release frequency.

## Reverse dependency checks

On 2026-11-10, we resolved the current CRAN reverse dependencies using tools::package_dependencies(reverse = TRUE, which = "most", recursive = "strong"): glydraw 0.9.0, glymotif 1.0.1, and glyparse 1.0.1.

We checked these CRAN packages against released glyrepr 1.2.0 and candidate 1.2.0.9000 using tools::check_packages_in_dir() with reverse dependencies enabled, and compared the results using tools::check_packages_in_dir_changes(). All six reverse-dependency checks completed with 0 errors, 0 warnings and 0 notes. There were no check-status regressions and no packages failed to be checked.
