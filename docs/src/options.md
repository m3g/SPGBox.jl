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


