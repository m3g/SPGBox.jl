# User guide

## Definition of the objective function and gradient functions

A function must be defined receiving as argument the current point as a vector: 

```julia-repl
julia> func(x) = x[1]^2 + x[2]^2

```

And the gradient must receive as arguments the vector of variables and a
vector which will be modified to contain the gradient at the current point:

```julia-repl
julia> function grad!(x,g)
         g[1] = 2*x[1]
         g[2] = 2*x[2]
       end

```
By Julia convention, to indicate that the gradient function modifies the vector `g`, we add
the `!` to its name, although this does not affect at all its behavior.

## Calling the solver, without bounds

The solver `spgbox!`, which modifies the input value of `x`, has a
minimal calling syntax of

```julia-repl
julia> x = rand(2)

julia> R = spgbox!(x,func,grad!)

```

The results will be returned to the data structure `R` of type
`SPGBoxResult`, and will be output as: 

```julia-repl
julia> R = spgbox!(x,func,grad!)

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 3.0
 Best solution found = [ 1.4210854715202004e-14, 2.0]
 Projected gradient norm = 1.4210854715202004e-13

 Number of iterations = 3
 Number of function evaluations = 5

```

## Calling the solver, with box bounds

Box bounds can be provided by the optional keywords `l` and `u` for
lower and upper bounds, respectively. These are vectors that will
delimit the bounds for each variable. For example, assuming the same
function and gradient functions defined in the example above, a lower
bound will be set for the second variable:

```julia-repl

julia> R = spgbox!(x,func,grad!,l=[-Inf,5])

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 12.0
 Best solution found = [ 0.0, 5.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 6

```

Upper bounds can be similarly set with `u=[+Inf,-5]`, for example.

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

### Constant parameters and new function and gradient methods 

The solver requires a function with a single argument, `x`, and a gradient
function with two arguments, `x` and `g`. If the function and gradient evalutions
require more parameters, use, for example: 

```julia-repl
julia> func(x,a,b,c) = a*x[1]^2 + (x[2]-b)^2 + c

julia> const a = 5. ; const b = 2. ; const c = 3. ;

julia> func(x) = func(x,a,b,c) 

```
To preserve performance it is fundamental to declare the parameters, in this
case `a`, `b`, and `c`, as constants (using `const`), to guarantee their
type-stability. This will allow the function specializations and compiler
optimizations that make Julia fast. 

The gradient function will be defined accordingly:

```julia-repl
julia> function grad!(x,g,a,b)
         g[1] = 2*a*x[1]
         g[2] = 2*(x[2]-b)
       end

julia> grad!(x,g) = grad!(x,g,a,b) 

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
julia> R = spgbox!(x, x -> func(x,a,b,c), (x,g) -> grad!(x,g,a,b))

```
where the second argument, `x -> func(x,a,b,c)` indicates that the objective
function is an anonymous function that, given `x`, returns `f(x,a,b,c)`. The gradient
is evaluated by an anonymous function that, given `(x,g)`, returns `grad!(x,g,a,b)`.  
This syntax also preserves performance and does not require the parameters to be declared
as constants. 

## Using automatic differentiation (AD)

Julia provides various packages for automatic differentiation, which can
be used in combination with `SPGBox`. Here, we illustrate the use of
[ForwardDiff](!http://www.juliadiff.org/ForwardDiff.jl/stable/). The
only two points that must be taken into consideration are: 1) The
AD function must modify an existing gradient vector and 2) use anonymous
closures to provide the gradient calculation function to the solver.

Here, a simple example, in which we use `ForwardDiff.gradient!` to
compute the derivative of a function which is the sum of squares of
the variables:

```julia-repl
julia> using SPGBox, ForwardDiff

julia> function func(x)
         f = 0.
         for i in 1:length(x)
           f += x[i]^2
         end
         f
       end

julia> x = rand(2)

julia> spgbox!(x,func,(x,g) -> ForwardDiff.gradient!(g,func,x) )

 SPGBOX RESULT:

 Convergence achieved.

 Final objective function value = 0.0
 Best solution found = [ 0.0, 0.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 3

```
