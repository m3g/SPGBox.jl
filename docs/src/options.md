# [Options](@id options)

Several keyword parameters can be used to adjust the convergence
criteria and to deal with the memory management of the execution.  

## Convergence criteria

Parameters exist to set the convergence threshold, maximum number of
iterations, maximum number of functional evaluations, and number of
possibly non-monotone steps.

These keywords provided to `spgbox!` with, for example:

``julia
julia> R = spgbox!(x,func,grad!,nitmax=1000)


where `nitmax`, in this case, is the maximum number of iterations.

The available keywords are:

| Keyword       | Type          |Meaning        | Default value  |
|:-------------:|:-------------:|:-------------:|:--------------:|
| `nitmax`      | `Integer` | Maximum number of iterations allowed. | `100` |
| `nfevalmax`   | `Integer` | Maximum number of function evaluations allowed. | `1000` |
| `eps`         | `Real` | Convergence criteria for the projected gradient norm. | `1e-5` |
| `m`           | `Integer` | Number of non-monotone search steps.  | `10` |

Additionally, a `iprint` option, with default value `0` is available.
Setting it to `1` or `2` will increase the detail of the iteration
output. 

## Memory management 

The SPGBox method requires four auxiliary vectors, three of them of
length equal to the number of variables (`g`, `xn` and `gn`), 
and a vector of length equal to
the number of non-monotone steps allowed, `fprev`, of dimension, `m`. 
Additionally, pre-allocation of the lower and upper bound arrays `l` an
`u` might be important, even if the problem is not bounded. 

For example, let us minimize the sum of squares of one million 
variables:

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

Without pre-allocating the auxiliary arrays:

```julia
julia> using BenchmarkTools

julia> @btime spgbox!(x,func,grad!)
  18.324 ms (13 allocations: 38.15 MiB)

```

Now we will pre-allocate all auxiliary arrays, and the lower and upper
bound arrays. Note that the `fprev` vector has dimension `10` in this
example, which is the default value for `m`. If `m` is modified, the
dimension of pre-allocated `fprev` must be set accordingly.

```julia
julia> g = zeros(n) ; xn = zeros(n); gn = zeros(n) ; fprev = zeros(10);

julia> l = [-Inf for _ in 1:n ]; u = [+Inf for _ in 1:n ];

```
And these arrays will be passed as arguments to the `spgbox!` function:

```julia
julia> @btime spgbox!(x,func,grad!,g=g,xn=xn,gn=gn,fprev=fprev,l=l,u=u)
  17.070 ms (4 allocations: 256 bytes)

```

While SPG is very memory efficient, pré-allocation of the arrays reduces
significanltly the use of memory, which might be important for multiple
executions of the the same code.



