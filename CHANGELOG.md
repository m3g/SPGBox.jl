SPGBox.jl Changelog
===================

[badge-feature]: https://img.shields.io/badge/Feature-green.svg
[badge-bugfix]: https://img.shields.io/badge/BUGFIX-yellow.svg
[badge-info]: https://img.shields.io/badge/Info-gray.svg

Version 0.7.4-DEV
-----------------
- ![FEATURE][badge-feature] Add the optional Optimization.jl interface through `SPGBoxOpt`, including box bounds and common solve options.
- ![BUGFIX][badge-bugfix] The iteration counter `nit` no longer counts the initial point: `nitmax = N` now allows `N` iterations (previously `N - 1`), and the reported `nit` is one less than before.
- ![INFO][badge-info] Document Optimization.jl usage and add a benchmark against other box-constrained solvers.

For older versions check the release log.
