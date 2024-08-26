""" 
    Vaux(x, fx; m=10)

A structure to store the auxiliary variables for the optimization algorithms.

It must be initialized with the current point `x`, and a value `fx` of the same
type of the output of the objective function. 

Optionally, the number of previous function values to store can be set with the
`m` keyword argument.

# Example

```jldoctest
julia> using SPGBox

julia> f(x) = x[1]^2 + x[2]^2 + 1;

julia> g!(g,x) = g .= 2 .* x;

julia> x0 = [5.0, 3.0];

julia> vaux = SPGBox.VAux(x0, 0.0; m=5);

julia> spgbox!(f, g!, x0; vaux=vaux)

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 1.0
 Sample of best point = Vector{Float64}[ 0.0, 0.0]
 Projected gradient norm = 0.0

 Number of iterations = 3
 Number of function evaluations = 3
```

"""
struct VAux{T,F}
    g::T
    xn::T
    gn::T
    fprev::Vector{F}
end
function VAux(x, fx; m=10)
    g = similar(x)
    xn = similar(x)
    gn = similar(x)
    fprev = [ zero(fx) for _ in 1:m ]
    return VAux{typeof(g),typeof(fx)}(g, xn, gn, fprev)
end