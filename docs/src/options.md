# [Options](@id Options)

Several keyword parameters can be used to adjust the convergence
criteria and to deal with the memory management of the execution.  

## Mutating or not mutating the initial point

Two solver functions are available: `spgbox!` and `spgbox`. The only difference is that `spgbox!` operates inplace in the input vector, while `spgbox` creates a copy of it and leaves the input variable unchanged. To perform completely non-allocating executions, use `spgbox!` and [preallocate the auxiliary vectors](https://m3g.github.io/SPGBox.jl/stable/options/#Memory-preallocation).

## Convergence criteria

Parameters exist to set the convergence threshold, maximum number of
iterations, maximum number of functional evaluations, and number of
possibly non-monotone steps.

These keywords provided to `spgbox!` with, for example:

```julia-repl
julia> R = spgbox!(f,g!,x,nitmax=1000)

```

where `nitmax`, in this case, is the maximum number of iterations.

The available keywords are:

| Keyword       | Type          |Meaning        | Default value  |
|:-------------:|:-------------:|:-------------:|:--------------:|
| `nitmax`      | `Integer` | Maximum number of iterations allowed. | `100` |
| `nfevalmax`   | `Integer` | Maximum number of function evaluations allowed. | `1000` |
| `eps`         | `eltype(x)` | Convergence criteria for the projected gradient norm. | `1e-5` |
| `m`           | `Integer` | Number of non-monotone search steps.  | `10` |

## Memory preallocation

The SPGBox method requires four auxiliary vectors, three of them of
length equal to the number of variables (`g`, `xn` and `gn`), 
and a vector of length equal to
the number of non-monotone steps allowed, `fprev`, of dimension, `m` 
(see [below](@ref fprev)). 

These vectors are allocated in the `SPGBox.VAux` data structure. For
preallocating the auxiliary vectors, initialize this data structure
before calling `spgbox!` and pass the data structure using the
`vaux` argument. In brief, given that the vector of variables is `x`, 
and the function is `f(x)`, do:

```julia-repl
julia> x = rand(1_000_000);

julia> auxvecs = SPGBox.VAux(x,f(x))

julia> R = spgbox!(f,g!,x,vaux=auxvecs)

```

For example, let us minimize the sum of squares of one million variables:

```julia-repl
julia> function f(x)
         f = 0.
         for i in eachindex(x)
           f += x[i]^2
         end
         f
       end
f (generic function with 1 methods)

julia> function g!(g,x)
         for i in eachindex(x)
           g[i] = 2*x[i]
         end
       end
g! (generic function with 1 methods)

julia> n = 1_000_000

julia> x = rand(n);

```

Without preallocating the auxiliary arrays:

```julia-repl
julia> using BenchmarkTools

julia> @btime spgbox!($f,$g!,$x)
  6.639 ms (10 allocations: 22.89 MiB)

```

Now we will preallocate all auxiliary arrays. 

```julia-repl
julia> auxvecs = SPGBox.VAux(x,f(x));

```
And these arrays will be passed as arguments to the `spgbox!` function:

```julia-repl
julia> @btime spgbox!($f,$g!,$x,vaux=$auxvecs)
  6.429 ms (0 allocations: 0 bytes)

```

While SPG is very memory efficient, pré-allocation of the arrays reduces
significanltly the use of memory, which might be important for multiple
executions of the the same code, for example in multi-start approach for
global-minimum search.

### [Size and preallocation of `fprev`](@id fprev)

The auxiliary vector `fprev` stores the information of the function value of the
last `m` function evaluations, which is the number of non-monotone
steps allowed. Thus, it is a vector of size `m` which is also
preallocated by `SPGBox.VAux`. By default, it is allocated to length
`10`, which is the default value for `m`. If `m` is modified by the
user and preallocation will be done, the new `m` value must be provided
to `SPGBox.VAux`, with:

```julia-repl
julia> auxvecs = SPGBox.VAux(x,f(x),m=5)

```

# Additional keywords available

Additional keywords available:

| Keyword       | Type          |Meaning        | Default value  |
|:-------------:|:-------------:|:-------------:|:--------------:|
| `iprint`      | `Integer` | Printing details (0, 1, or 2) | `0` |
| `project_x0`  | `Bool` | Projects, or not, the initial point on the bounds. | `true` |




