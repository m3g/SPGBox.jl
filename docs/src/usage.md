# Usage

## Definition of the objective function and gradient functions

A function must be defined receiving as arguments the current point as a
vector. If the function requires additional parameters, these can
provided by defining an additional parameter-free method. For example: 

```julia
julia> func(x,a,b,c) = a*x[1]^2 + (x[2]-b)^2 + c

julia> a = 5. ; b = 2. ; c = 3. ;

julia> func(x) = func(x,a,b,c) 

```

A function to compute the gradient must also be provided. This function
will assume that the gradient vector is already allocated and will
modify it. The parameters of the gradient function must be at least,
therefore, the current point `x` and the gradient vector `g`:

```julia
julia> function grad!(x,g,a,b,c)
         g[1] = 2*a*x[1]
         g[2] = 2*(x[2]-b)
       end

julia> grad!(x,g) = grad!(x,g,a,b,c) 

```

The function method which receive only the current point `x`, and the
gradient method which receives only `x` and the gradient vector `g` are
the ones actually invoked by the optimizer. By convention, the gradient
function name includes an `!` to indicate that it modifies its
arguments.

## Calling the optimizer, without bounds

The optimizer `spgbox!`, which modifies the input value of `x`, has a
minimal calling syntax of

```julia
julia> x = rand(2)

julia> R = spgbox!(x,func,grad!)

```

The results will be returned to the data structure `R` of type
`SPGBoxResult`, and will be output as: 

```julia
julia> R = spgbox!(x,func,grad!)

julia> R = spgbox!(x,func,grad!)

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 3.0
 Best solution found = [ 1.4210854715202004e-14, 2.0]
 Projected gradient norm = 1.4210854715202004e-13

 Number of iterations = 3
 Number of function evaluations = 5

```

## Calling the optimizer, with box bounds

Box bounds can be provided by the optional keywords `l` and `u` for
lower and upper bounds, respectively. These are vectors that will
delimit the bounds for each variable. For example, assuming the same
function and gradient functions defined in the example above, a lower
bound will be set for the second variable:

```julia

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

```julia
struct SPGBoxResult
  x :: Vector{Float64}
  f :: Float64
  gnorm :: Float64
  nit :: Int64
  nfeval :: Int64
  ierr :: Int64
end
```

The data can be accesses as usual, using, for example:

```julia
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

|                  |                 |
|:-----------------|:----------------|
| `ierr=0`           | Success: convergence achieved. |
| `ierr=1`           | Maximum number of iterations achieved. |
| `ierr=2`           | Maximum number of function evaluations achieved.  |

The convergence criteria can be adjusted using optional keywords, as
described in the [Options](@ref Options) section.






