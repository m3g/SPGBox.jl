# User guide

## Definition of the objective function and gradient functions

A function must be defined receiving as argument the current point as a vector: 

```julia-repl
julia> f(x) = x[1]^2 + (x[2]-2)^2
```

And the gradient must receive as arguments the vector of variables and a
vector which will be modified to contain the gradient at the current point:

```julia-repl
julia> function g!(g,x)
         g[1] = 2*x[1]
         g[2] = 2*(x[2]-2)
       end
```
By Julia convention, to indicate that the gradient function modifies the vector `g`, we add
the `!` to its name, although this does not affect at all its behavior.

## Calling the solver, without bounds

The solver function is `spgbox!`, which mutates the input value of `x` (with the best solution found in it at the output). Alternatively, use `spgbox` to call the solver without modifying the values of the input vector (although a copy of if this vector will be produced). 

The solver calls have a minimal calling syntax of
```julia-repl
julia> x = rand(2);

julia> R = spgbox!(f,g!,x)
```

The results will be returned to the data structure `R` of type
`SPGBoxResult`, and will be output as: 

```julia-repl
julia> R = spgbox!(f,g!,x)

 SPGBOX RESULT:

 Convergence achieved.

 Final objective function value = 0.0
 Best solution found = [ 0.0, 2.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 3
```

## Calling the solver, with box bounds

Box bounds can be provided by the optional keywords `lower` and `upper` for
lower and upper bounds, respectively. These are vectors that will
delimit the bounds for each variable. For example, assuming the same
function and gradient functions defined in the example above, a lower
bound will be set for the second variable:

```julia-repl
julia> x = rand(2);

julia> R = spgbox!(f,g!,x,lower=[-Inf,5])

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 9.0
 Best solution found = [ 0.0, 5.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 3
```

Upper bounds can be similarly set with `upper=[+Inf,-5]`, for example.

Note, the bounds can also be provided as non-keyword parameters, with:
```
julia> lower = [-Inf,5]; upper = [+Inf, -2];

julia> R = spgbox!(f,g!,x,lower,upper)
```

## Result data structure and possible outcomes

The minimization can result in a successful convergence, on in exceeded
numbers of iterations or functional evaluations. These outcomes are
explicit in the output printed (second line), and stored in the result
structure, which contains the following data: 

```julia-repl
struct SPGBoxResult
  x :: Vector{Float64}
  f :: Float64
  gnorm :: Float64
  nit :: Int64
  nfeval :: Int64
  ierr :: Int64
end
```

The data can be accessed as usual, using, for example:

```julia-repl
julia> R.f
12.0
```

The data structure contains:

| Variable name    | Meaning         |
|:-----------------|:----------------|
| `x`              | Best point found (solution if `ierr=0`) |
| `f`              | Best function value found. |
| `gnorm`          | Norm of gradient projected on the constraints. |
| `nit`            | Number of iterations performed. | 
| `nfeval`         | Number of function evaluations. |
| `ierr`           | Exit status.  |

The possible outcomes of `ierr` are:

| Outcome          |  Meaning                 |
|:-----------------|:----------------|
| `ierr=0`           | Success: convergence achieved. |
| `ierr=1`           | Maximum number of iterations achieved. |
| `ierr=2`           | Maximum number of function evaluations achieved.  |

The convergence criteria can be adjusted using optional keywords, as
described in the [Options](@ref Options) section.

## Data-dependent function evaluation

If the function requires additional parameters, two strategies are
possible while preserving performance: 1) Declare the parameters as constants
and define an extra method, or 2) Pass the function as an anonymous closure. 

## Input data types

The implementation of `SPGBox` is generic for input types of any abstract matrix or vector. 
Propagation of units, or numbers with arbitrary precision is supported. As minimal example illustrating 
these features, let us minimize a function where the input variables are given as a matrix of `BigFloat`s,
with units:

```julia-repl
julia> function f(x)
           f = zero(typeof(x[begin]^2))
           for i in eachindex(x)
               f += (x[i] - i*oneunit(eltype(x)))^2
           end
           return f
       end
f (generic function with 1 method)

julia> function g!(g,x)
           g .= zero(eltype(x))
           for i in eachindex(x)
               g[i] = 2*(x[i] - i*oneunit(eltype(x)))
           end
           return g
       end
g! (generic function with 1 method)

julia> x = rand(BigFloat,2,2)u"nm"
2×2 Matrix{Quantity{BigFloat, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}:
 0.128083 nm  0.817173 nm
 0.139545 nm  0.391047 nm

julia> spgbox(f,g!,x)

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 0.0 nm^2
 Sample of best point = Matrix{Quantity{BigFloat, 𝐋, Unitful.FreeUnits{(nm,), 𝐋, nothing}}}[ 1.0 nm, 2.0 nm, 3.0 nm, 4.0 nm]
 Projected gradient norm = 0.0 nm

 Number of iterations = 2
 Number of function evaluations = 3

```

### Constant parameters and new function and gradient methods 

The solver requires a function with a single argument, `x`, and a gradient
function with two arguments, `x` and `g`. If the function and gradient evalutions
require more parameters, use, for example: 

```julia-repl
julia> f(x,a,b,c) = a*x[1]^2 + (x[2]-b)^2 + c

julia> const a = 5. ; const b = 2. ; const c = 3. ;

julia> f(x) = f(x,a,b,c) 
```
To preserve performance it is fundamental to declare the parameters, in this
case `a`, `b`, and `c`, as constants (using `const`), to guarantee their
type-stability. This will allow the function specializations and compiler
optimizations that make Julia fast. 

The gradient function will be defined accordingly:

```julia-repl
julia> function g!(g,x,a,b)
         g[1] = 2*a*x[1]
         g[2] = 2*(x[2]-b)
       end

julia> g!(g,x) = g!(g,x,a,b) 
```

The function method which receives only the current point `x`, and the
gradient method which receives only `x` and the gradient vector `g` are
the ones actually invoked by the solver.

### Using anonymous closures 

An anonymous closure is a function with a special syntax of the form 

```julia
x -> f(x)
```
which should be read as "given `x`, return `f(x)`". These anonymous functions can
be provided directly as arguments to the solver, while providing an interface for 
using external parameters. Considering the same function and gradient functions
above, one uses anonymous functions  directly as arguments in the solver call:

```julia-repl
julia> R = spgbox!(x -> f(x,a,b,c), (g,x) -> g!(g,x,a,b), x)
```
where the first argument, `x -> f(x,a,b,c)` indicates that the objective
function is an anonymous function that, given `x`, returns `f(x,a,b,c)`. The gradient
is evaluated by an anonymous function that, given `(g,x)`, returns `g!(g,x,a,b)`.  
This syntax also preserves performance and does not require the parameters to be declared
as constants. 

## Using automatic differentiation (AD)

Julia provides various packages for automatic differentiation, which can
be used in combination with `SPGBox`. Here, we illustrate the use of
[ReverseDiff](http://www.juliadiff.org/ReverseDiff.jl/). The
only two points that must be taken into consideration are: 1) The
AD function must modify an existing gradient vector and 2) use anonymous
closures to provide the gradient calculation function to the solver.

Here, a simple example, in which we use `ReverseDiff.gradient!` to
compute the derivative of a function which is the sum of squares of
the variables:

```julia-repl
julia> using SPGBox, ReverseDiff

julia> function f(x)
         f = 0.
         for i in eachindex(x)
           f += x[i]^2
         end
         f
       end

julia> x = rand(2)

julia> spgbox!(f, (g,x) -> ReverseDiff.gradient!(g,f,x), x, lower=[-Inf,2.])

 SPGBOX RESULT:

 Convergence achieved.

 Final objective function value = 4.0
 Best solution found = [ 0.0, 2.0]
 Projected gradient norm = 0.0

 Number of iterations = 0
 Number of function evaluations = 1
```
