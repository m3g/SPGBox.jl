# SPGBox

Implements the Spectral Projected Gradient Method for minimization 
with box constraints. 

Algorithm reference:

E. G. Birgin, J. M. Martínez and M. Raydan, "Nonmonotone spectral
projected gradient methods on convex sets", SIAM Journal on Optimization
10, pp. 1196-1211, 2000. 
[LINK](http://www.ime.usp.br/~egbirgin/publications/bmr.pdf)

## Documentation

The complete documentation can be found at: 
https://m3g.github.io/SPGBox.jl/stable/

## How to install

```julia-repl
julia> using Pkg

julia> Pkg.add("SPGBox")
```

or, more concisely,

```julia-repl
julia> ] add SPGBox
```

## Quick usage example:

Define the function to compute the objective function and the gradient,
for example with:

```julia-repl
julia> f(x) = x[1]^2 + x[2]^2

julia> function g!(x,g)
         g[1] = 2*x[1]
         g[2] = 2*x[2]
       end
```

And the objective function can be minimized with optional box bounds.
Here, with a lower bound of `2` for the first variable:

```julia-repl
julia> x = 2 .+ rand(2)

julia> spgbox!(f,g!,x,lower=[2.,-Inf])

 SPGBOX RESULT:

 Convergence achieved.

 Final objective function value = 4.0
 Best solution found = [ 2.0, 0.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 3

```

