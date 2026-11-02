## Package maintenance

This is the first release maintained by Laszlo Rado <lazlorado@gmail.com>
following the transfer of package maintenance from Stuart Lacy
<stuart.lacy@york.ac.uk>. Stuart Lacy remains listed as an author.

Stuart Lacy concluded his maintenance of the package, and Laszlo Rado
volunteered to take over its continued development and CRAN maintenance.
Stuart Lacy confirmed by email that he supports continued development of the
fork and its submission to CRAN under the `rprev` package name. The original
confirmation can be provided to CRAN on request.

The package was archived on CRAN on 2026-03-18 because calls to deprecated
`dplyr` APIs had not been updated. Version 1.1.0 replaces those calls with
their current equivalents and restores compatibility with current dependency
versions.

## Release summary

Version 1.1.0 adds prevalence estimation for one or more index dates in a single
call. Within each bootstrap iteration, all requested dates use the same
simulated incident population and individual-level survival threshold, producing
temporally consistent survival states. Single-index seeded simulation behaviour
is retained.

The release also updates per-date outputs, S3 methods, documentation, dependency
compatibility, and unit test coverage.

## Test environments

* local Windows 11 x86_64 (build 26100), R 4.6.2 (2026-10-29 ucrt),
  `R CMD check --as-cran --no-manual`
* win-builder Windows Server 2022 x86_64, R-devel (2026-10-31 r90820 ucrt)
* R-hub Ubuntu 26.04.1 LTS x86_64, R-devel
* R-hub macOS Sequoia 15.7.9 x86_64, R-devel

## R CMD check results

* local Windows: 0 errors | 0 warnings | 1 note
* win-builder R-devel: 0 errors | 0 warnings | 1 note
* R-hub Linux R-devel: 0 errors | 0 warnings | 0 notes
* R-hub macOS R-devel: 0 errors | 0 warnings | 0 notes

The NOTE identifies this as a new submission and reports that the previous CRAN
release was archived. The archive reason and maintenance transfer are addressed
above.

The local check used `--no-manual` because the available TinyTeX 2022
installation does not contain `makeindex` and cannot install 2026 TeX Live
components across releases. The PDF and HTML manuals were successfully checked
by win-builder R-devel.

The complete `testthat` suite passes with 230 tests and no failures, warnings,
or skips on local Windows, R-hub Linux, and R-hub macOS.

## Downstream dependencies

The current CRAN package index lists no reverse dependencies for rprev.
