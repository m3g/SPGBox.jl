const SPGBOXOPT_OPTIONS = (:eps, :nfevalmax, :m, :project_x0, :step_nc, :iprint)

"""
    SPGBoxOpt(; eps, nfevalmax, m, project_x0, step_nc, iprint)

Algorithm to use SPGBox through the [Optimization.jl](https://docs.sciml.ai/Optimization/stable/) interface.

The solver is available once `Optimization` (or `OptimizationBase` and `SciMLBase`) is loaded:

```julia
using SPGBox, Optimization, ForwardDiff
f(x, p) = (x[1] - p[1])^2 + (x[2] - 1)^2
optf = OptimizationFunction(f, AutoForwardDiff())
prob = OptimizationProblem(optf, [10.0, 18.0], [0.0]; lb = [2.0, -Inf], ub = [Inf, Inf])
sol = solve(prob, SPGBoxOpt())
```

The keyword arguments are the options of [`spgbox!`](@ref) of the same name, and are the 
only ones that can be set in the constructor. The common `solve` options are interpreted as:

- `maxiters`: sets `nitmax`. If set, the limit on the number of function evaluations is
  `nfevalmax` if provided to `SPGBoxOpt`, and unlimited otherwise.
  If not set, the defaults of `spgbox!` are used.
- `maxtime`: maximum time (in seconds), checked at the end of every iteration.
- `callback`: called as `callback(state, f)` after every iteration. The `state.original` field is the 
  current `SPGBoxResult`. If the callback returns `true` the optimization stops.
- `abstol` and `reltol` are not used. Use the `eps` option (tolerance on the projected gradient norm).

The gradient is computed by the automatic differentiation backend of the `OptimizationFunction`,
or provided with its `grad` option.
"""
struct SPGBoxOpt{K<:NamedTuple}
    kwargs::K
end

function SPGBoxOpt(; kwargs...)
    for k in keys(kwargs)
        k in SPGBOXOPT_OPTIONS || throw(ArgumentError(
            "Unknown option `$k` for SPGBoxOpt. Valid options: $(join(SPGBOXOPT_OPTIONS, ", "))."
        ))
    end
    return SPGBoxOpt(NamedTuple(kwargs))
end

Base.show(io::IO, opt::SPGBoxOpt) =
    print(io, "SPGBoxOpt(", join(("$k = $(repr(v))" for (k, v) in pairs(opt.kwargs)), ", "), ")")
