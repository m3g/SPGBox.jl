# Optimization.jl benchmark

This benchmark compares `SPGBox` (through its `SPGBoxOpt` [Optimization.jl](https://docs.sciml.ai/Optimization/stable/)
interface) with the other gradient-based, box-constrained *local* solvers available
through Optimization.jl:

- `OptimizationLBFGSB.LBFGSB()`: the classical L-BFGS-B algorithm.
- `Optim.Fminbox(Optim.LBFGS())`: L-BFGS with box constraints enforced by an
  interior/barrier approach.
- `NLopt.LD_LBFGS()`: NLopt's L-BFGS, which supports bounds directly.

These are the solvers most comparable to SPGBox: deterministic, gradient-based,
local methods that natively support box constraints (unlike the global/stochastic
solvers compared in the
[SciML black-box benchmark](https://docs.sciml.ai/SciMLBenchmarksOutput/dev/GlobalOptimization/blackbox_global_optimizers/),
which target non-smooth or multimodal problems).

SPGBox implements a first-order method (spectral projected gradient), so it is
expected to require more function/gradient evaluations than the quasi-Newton
solvers above, in exchange for a much lower cost per iteration and no need to
store or update an approximate Hessian. The benchmark illustrates this
trade-off on three scalable test problems (an extended Rosenbrock function and
a diagonal quadratic, both with bounds that force part of the solution onto
the boundary, and the Powell singular function, which is ill-conditioned at
the solution) at dimensions 20, 100 and 1000.

## Running

```bash
cd benchmark/optimization
julia --project=. -e 'import Pkg; Pkg.develop(path="../.."); Pkg.instantiate()'
julia --project=. optimization_benchmark.jl
```

The `Pkg.develop` step is only needed to benchmark an unreleased version of
`SPGBox` (such as a branch under development); otherwise the registered
version is used.

This prints a results table to the terminal and writes it to `results.csv`,
with the final objective value, the number of function/gradient evaluations
(counted independently for every solver, since not all of them report this
through `sol.stats`), and the median solve time, for every
(problem, dimension, solver) combination.
