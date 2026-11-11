```@meta
CollapsedDocStrings = true
```

# [Optimization.jl interface](@id Optimization)

`SPGBox` can be called through the common
[Optimization.jl](https://docs.sciml.ai/Optimization/stable/) interface, using
the `SPGBoxOpt` algorithm. This allows `SPGBox` to be used interchangeably
with other Optimization.jl solvers, and to be compared with them on the same
problems.

The interface is provided by a package extension, loaded automatically when
`SPGBox` and `Optimization` (or `OptimizationBase` and `SciMLBase`) are both
loaded. No extra dependency is added to `SPGBox` for users who do not need
this interface.

## Basic usage

The objective function must be of the Optimization.jl form `f(x,p)`, where
`p` are optional fixed parameters, and the gradient is computed with an
automatic-differentiation backend (here, `ForwardDiff`):

```julia-repl
julia> using SPGBox, Optimization, ForwardDiff

julia> f(x, p) = (x[1] - p[1])^2 + (x[2] - 1)^2
f (generic function with 1 method)

julia> optf = OptimizationFunction(f, AutoForwardDiff())
OptimizationFunction{...}

julia> prob = OptimizationProblem(optf, [10.0, 18.0], [1.0]; lb = [2.0, -Inf], ub = [Inf, Inf])
OptimizationProblem{...}

julia> sol = solve(prob, SPGBoxOpt())
retcode: Success
u: [2.0, 1.0]
Final objective value:     1.0
```

A gradient function can be provided instead of an AD backend, with the
`grad` keyword of `OptimizationFunction` (see the
[Optimization.jl documentation](https://docs.sciml.ai/Optimization/stable/API/optimization_function/)).

## `SPGBoxOpt` options

```@docs
SPGBoxOpt
```

## Mapping of the common `solve` options

| `solve` keyword | Effect |
|:----------------|:-------|
| `maxiters`      | Sets `nitmax`. If given, the evaluation limit is the `nfevalmax` passed to `SPGBoxOpt` (or unlimited, if none was passed). If not given, the `spgbox!` defaults are used. |
| `maxtime`       | Checked at the end of every iteration (in seconds). |
| `callback`      | Called as `callback(state, f)` after every iteration. `state.original` holds the current `SPGBoxResult`. Returning `true` stops the optimization. |
| `abstol`, `reltol` | Not used. Use the `eps` option of `SPGBoxOpt` (tolerance on the projected gradient norm) instead. |

The `lb` and `ub` bounds of the `OptimizationProblem` are passed directly to
`spgbox!`.

## Benchmark against other box-constrained solvers

A script comparing `SPGBoxOpt` against other gradient-based, box-constrained
Optimization.jl solvers (L-BFGS-B, `Fminbox(LBFGS())` and `NLopt.LD_LBFGS()`)
is available in
[`benchmark/optimization/`](https://github.com/m3g/SPGBox.jl/tree/main/benchmark/optimization)
in the repository.
