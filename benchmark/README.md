# Benchmark

This code may be used to benchmark the current version of `SPGBox.jl`. Just move
here and run
```bash
julia --project=. benchmark.jl
```
It will use `SPGBox.jl` to solve a battery of tests from CUTER similar to the
tests used in the original SPG paper. After a while it will generate a file
named `bench_results.csv` that has performance metrics for solving each test
problem like number of iterations (NITER), number of funcion evaluations
(FEVAL), elapsed time (TIME). If these metrics are negative is because the
solved failed.

It also has a simple tool that will generate a performace profile from two files
named `base_bench.csv` and `bench_results.csv`. You can do this running
```bash
julia --project=. genpprof.jl FEVAL
```
This will generate the profile that compare the runs with respect to the
number of function evaluations in a file named `profile.png`. You can change
FEVAL for NITER ou TIME.
