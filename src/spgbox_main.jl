#
# Algorithm of:
#
# NONMONOTONE SPECTRAL PROJECTED GRADIENT METHODS ON CONVEX SETS
# ERNESTO G. BIRGIN, JOSÉ MARIO MARTÍNEZ, AND MARCOS RAYDAN
# SIAM J. O. PTIM. Vol. 10, No. 4, pp. 1196-1211
#
# Implemented by J. M. Martínez (IMECC - UNICAMP)
# Iinitally translated and adapted to Julia by L. Martínez (IQ-UNICAMP)
#
"""
```
spgbox!(f, g!, x::AbstractVecOrMat; lower=..., upper=..., options...)`
```

```
spgbox!(f, g!, lower::AbstractVeOrMat, upper::AbstractVecOrMat, x::AbstractVecOrMat; options...)`
```

Minimizes function `f` starting from initial point `x`, given the function to compute the gradient, `g!`. `f` must be of the form `f(x)`, and `g!` of the form `g!(g,x)`, where `g` is the gradient vector to be modified. It modifies the `x` vector, which will contain the best solution found (see `spgbox` for a non-mutating alternative). 

Optional lower and upper box bounds can be provided using optional arguments `lower` and `upper`, which can be provided as the fourth and fifth arguments or with keyword parameters.

Returns a structure of type `SPGBoxResult`, containing the best solution found in `x` and the final objective function in `f`.

# Examples
```jldocstest
julia> f(x) = x[1]^2 + x[2]^2

julia> function g!(g,x)
           g[1] = 2*x[1]
           g[2] = 2*x[2]
       end
```
## Without bounds

```jldocstest
julia> x = rand(2)

julia> spgbox!(f,g!,x)

SPGBOX RESULT: 

Convergence achieved. 

Final objective function value = 0.0
Sample of best point = Vector{Float64}[ 0.0, 0.0]
Projected gradient norm = 0.0

Number of iterations = 3
Number of function evaluations = 3

```

## With bounds

```jldocstest
julia> x = 2 .+ rand(2)

julia> spgbox!(f,g!,x,lower=[2.,-Inf])

SPGBOX RESULT: 

Convergence achieved. 

Final objective function value = 4.0
Sample of best point = Vector{Float64}[ 2.0, 0.0]
Projected gradient norm = 0.0

Number of iterations = 3
Number of function evaluations = 3

```
"""
function spgbox!(f::F, g!::G, x::AbstractVecOrMat{T}; kargs...) where {F<:Function,G<:Function,T} 
    @show f(x)
    spgbox!(
        (g,x) -> begin 
            g!(g,x)
            return f(x)
        end,
        x, f; kargs...
    )
end

#
# Call with a single function to compute the function and the gradient
#
function spgbox!(
    fg!::Function,
    x::AbstractVecOrMat{T},
    func_only = nothing;
    lower::Union{Nothing,AbstractVecOrMat{T}} = nothing,
    upper::Union{Nothing,AbstractVecOrMat{T}} = nothing,
    eps = 1.e-5 * oneunit(T),
    nitmax::Int = 100,
    nfevalmax::Int = 1000,
    m::Int = 10,
    vaux::VAux = VAux(x, (isnothing(func_only) ? fg!(g,x) : func_only(x)), m=m),
    iprint::Int = 0,
    project_x0::Bool = true
) where {T}

    # Number of variables
    n = length(x)

    # Auxiliary arrays (associate names and check dimensions)
    g = vaux.g
    xn = vaux.xn
    gn = vaux.gn
    fprev = vaux.fprev
    @assert length(g) == n "Auxiliar gradient vector `g` must be of the same length as `x`"
    @assert length(xn) == n "Auxiliar vector `xn` must be of the same length as `x`"
    @assert length(gn) == n "Auxiliar vector `gn` must be of the same length as `x`"
    @assert length(fprev) == m "Auxiliar vector `fprev` must be of length `m`"

    # Check if bounds are defined, project or not the initial point on them
    if !isnothing(lower)
        @assert length(lower) == n "Lower bound vector `lower` must be of the same length than x, got: $(length(lower))"
        if project_x0
            @. x = max(x, lower)
        else
            for i in eachindex(x)
                if x[i] < lower[i]
                    error(
                        " Initial value of variable $i smaller than lower bound, and `project_x0` is set to `false`. ",
                    )
                end
            end
        end
    end
    if !isnothing(upper)
        @assert length(upper) == n "Upper bound vector `upper` must be of the same length than `x`, got: $(length(upper))"
        if project_x0
            @. x = min(x, upper)
        else
            for i in eachindex(x)
                if x[i] > lower[i]
                    error(
                        " Initial value of variable $i greater than upper bound, and `project_x0` is set to `false`. ",
                    )
                end
            end
        end
    end

    # Iteration counter
    nit = 1

    # Objective function and gradient at initial point
    nfeval = 1
    fcurrent = fg!(g,x)
    gnorm = pr_gradnorm(g, x, lower, upper)
    gnorm <= eps && return SPGBoxResult(x, fcurrent, gnorm, nit, nfeval, 0)

    # Initialize array of previous function values
    tspg = one(T)
    for i in eachindex(fprev)
        fprev[i] = fcurrent
    end

    while nit < nitmax

        if iprint > 0
            println("----------------------------------------------------------- ")
            println(" Iteration: ", nit)
            println(" x = ", x[begin], " ... ", x[end])
            println(" Objective function value = ", fcurrent)
            println(" ")
            println(" Norm of the projected gradient = ", gnorm)
            println(" Number of function evaluations = ", nfeval)
        end

        t = tspg
        fref = maximum(fprev)
        if iprint > 2
            println(" fref = ", fref)
            println(" fprev = ", fprev)
            println(" t = ", t)
        end

        # Evaluate trial point
        compute_xn!(xn,x,t,g,lower,upper)
        nfeval = nfeval + 1
        nfeval > nfevalmax && return SPGBoxResult(x, fcurrent, gnorm, nit, nfeval, 2)
        fn = fg!(gn, xn)

        while (fn > fref)
            if iprint > 2
                println(" xn = ", xn[begin], " ... ", xn[end])
                println(" f[end] = ", fn, " fref = ", fref)
            end
            compute_xn!(xn,x,t,g,lower,upper)
            nfeval = nfeval + 1
            nfeval > nfevalmax && return SPGBoxResult(x, fcurrent, gnorm, nit, nfeval, 2)
            if !isnothing(func_only)
                fn = func_only(xn)
            else
                fn = fg!(gn,xn)
            end
            # Reduce step
            t = t / 2
        end

        # Trial point accepted
        num = zero(T)
        den = zero(T)
        for i in eachindex(x)
            num = num + (xn[i] - x[i])^2 / oneunit(T)
            den = den + (xn[i] - x[i]) * (gn[i] - g[i]) / oneunit(T)
        end
        if den <= zero(T)
            tspg = T(100)
        else
            tspg = min(1000 * one(T), one(T) * num / den)
        end
        fcurrent = fn
        for i in eachindex(x)
            x[i] = xn[i]
            g[i] = gn[i]
        end
        for i = firstindex(fprev):lastindex(fprev)-1
            fprev[i] = fprev[i+1]
        end
        fprev[end] = fcurrent
        nit = nit + 1

        # Compute projected gradient norm
        gnorm = pr_gradnorm(g, x, lower, upper)
        gnorm <= eps && return SPGBoxResult(x, fcurrent, gnorm, nit, nfeval, 0)
    end

    # Maximum number of iterations achieved
    return SPGBoxResult(x, fcurrent, gnorm, nit, nfeval, 1)
end
#
# Call with lower and upper as positional arguments
#
spgbox!(f::F, g!::G, lower, upper, x; kargs...) where {F,G} = spgbox!(f, g!, x, lower=lower, upper=upper, kargs...)
spgbox!(fg!::F, lower, upper, x; kargs...) where {F,G} = spgbox!(fg!, x, lower=lower, upper=upper, kargs...)

"""

```
spgbox(f, g!, x::AbstractVecOrMat; lower=..., upper=..., options...)`
```

```
spgbox(f, g!, lower::AbstractVecOrMat, upper::AbstractVecOrMat, x::AbstractVecOrMat; options...)`
```

Minimizes function `f` starting from initial point `x`, given the function to compute the gradient, `g!`. `f` must be of the form `f(x)`, and `g!` of the form `g!(g,x)`, where `g` is the gradient vector to be modified. This call *does not* mutate the `x` vector, instead it will create a (deep)copy of it (see `spgbox!` to the inplace alternative). 

Optional lower and upper box bounds can be provided using optional arguments `lower` and `upper`, which can be provided as the fourth and fifth arguments or with keyword parameters.

Returns a structure of type `SPGBoxResult`, containing the best solution found in `x` and the final objective function in `f`.


"""
function spgbox(f::F, g!::G, x::AbstractVecOrMat{T}; kargs...) where {F<:Function,G<:Function,T}
    x0 = deepcopy(x)
    return spgbox!(f, g!, x0; kargs...)
end
# With a single function to compute the function and the gradient
function spgbox(fg::FG, x::AbstractVecOrMat{T}; kargs...) where {FG<:Function,T}
    x0 = deepcopy(x)
    return spgbox!(fg, x0; kargs...)
end
#
# call with lower and upper as positional arguments
#
function spgbox(f::F, g!::G, lower, upper, x::AbstractVecOrMat{T}; kargs...) where {F<:Function,G<:Function,T}
    x0 = deepcopy(x)
    return spgbox!(f, g!, lower, upper, x0; kargs...)
end
# With a single function to compute the function and the gradient
function spgbox(fg::FG, lower, upper, x::AbstractVecOrMat{T}; kargs...) where {FG<:Function,T}
    x0 = deepcopy(x)
    return spgbox!(fg, lower, upper, x0; kargs...)
end
