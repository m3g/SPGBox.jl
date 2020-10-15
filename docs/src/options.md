# [Options](@id Options)

Several keyword parameters can be used to adjust the convergence
criteria and to deal with the memory management of the execution.  

## Convergence criteria

Parameters exist to set the convergence threshold, maximum number of
iterations, maximum number of functional evaluations, and number of
possibly non-monotone steps.

These keywords provided to `spgbox!` with, for example:

```julia
julia> R = spgbox!(x,func,grad!,nitmax=1000)

```

where `nitmax`, in this case, is the maximum number of iterations.

The available keywords are:

| Keyword       | Type          |Meaning        | Default value  |
|:-------------:|:-------------:|:-------------:|:--------------:|
| `nitmax`      | `Integer` | Maximum number of iterations allowed. | `100` |
| `nfevalmax`   | `Integer` | Maximum number of function evaluations allowed. | `1000` |
| `eps`         | `Real` | Convergence criteria for the projected gradient norm. | `1e-5` |
| `m`           | `Integer` | Number of non-monotone search steps.  | `10` |

## Memory management 

The SPGBox method requires four auxiliary vectors, three of them of
length equal to the number of variables (`g`, `xn` and `gn`), 
and a vector of length equal to
the number of non-monotone steps allowed, `fprev`, of dimension, `m` 
(see [below](@ref fprev)). 

These vectors are allocated in the `SPGBox.Aux` data structure. For
preallocating the auxiliary vectors, initialize this data structure
before calling `spgbox!` and pass the data structure using the
`aux` argument. In brief, do, for 1 million variables:

```julia
julia> n = 1_000_000

julia> auxvecs = SPGBox.Aux(n)

julia> R = spgbox!(x,func,grad!,aux=auxvecs)

```

For example, let us minimize the sum of squares of one million variables:

```julia
julia> function func(x)
         f = 0.
         for i in 1:length(x)
           f += x[i]^2
         end
         f
       end
func (generic function with 2 methods)

julia> function grad!(x,g)
         for i in 1:length(x)
           g[i] = 2*x[i]
         end
       end
grad! (generic function with 2 methods)

julia> n = 1_000_000

julia> x = rand(n);

```

Without preallocating the auxiliary arrays:

```julia
julia> using BenchmarkTools

julia> @btime spgbox!($x,func,grad!)
  6.639 ms (10 allocations: 22.89 MiB)

```

Now we will preallocate all auxiliary arrays. 

```julia
julia> auxvecs = SPGBox.Aux(n);

```
And these arrays will be passed as arguments to the `spgbox!` function:

```julia
julia> @btime spgbox!($x,func,grad!,aux=$auxvecs)
  6.429 ms (0 allocations: 0 bytes)

```

While SPG is very memory efficient, pré-allocation of the arrays reduces
significanltly the use of memory, which might be important for multiple
executions of the the same code.

### [Size and preallocation of `fprev`](@id fprev)

The auxiliary vector `fprev` stores the information of the function value of the
last `m` function evaluations, which is the number of non-monotone
steps allowed. Thus, it is a vector of size `m` which is also
preallocated by `SPGBox.Aux`. By default, it is allocated to length
`10`, which is the default value for `m`. If `m` is modified by the
user and preallocation will be done, the new `m` value must be provided
to `SPGBox.Aux`, with:

```julia
julia> auxvecs = SPGBox.Aux(n,m)

```

# Additional keywords available

Additional keywords available:

| Keyword       | Type          |Meaning        | Default value  |
|:-------------:|:-------------:|:-------------:|:--------------:|
| `iprint`      | `Integer` | Printing details (0, 1, or 2) | `0` |
| `project_x0`  | `Bool` | Projects, or not, the initial point on the bounds. | `true` |




