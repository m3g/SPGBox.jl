"""
    SPGBoxResult

Data structure that contains the results of the `spgbox` optimization method. 

    `x`: Array containing final point (solution or best point found).
 
    `f`: Final objective function value.

    `gnorm`: Norm of the projected gradient on the bounds.

    `nit`: Number of iterations performed. 

    `nfeval`: Number of function evaluations.

    `ierr`: Status of execution: 

         0. Converged successfully to solution. 

         1. Maximum number of iterations achieved.

         2. Maximum number of function evaluations achieved.

    `return_from_callback`: Boolean indicating if the algorithm returned from the callback function.

"""
struct SPGBoxResult{T,F,G}
    x::T
    f::F
    gnorm::G
    nit::Int64
    nfeval::Int64
    ierr::Int64
    return_from_callback::Bool
end

function Base.show(io::IO, R::SPGBoxResult)
    println(io, "")
    println(io, " SPGBOX RESULT: ")
    println(io, "")
    if R.return_from_callback
        println(io, " Returned from callback function. ")
    else
        R.ierr == 0 && println(io, " Convergence achieved. ")
        R.ierr == 1 && println(io, " Maximum number of iterations (nitmax) reached.")
        R.ierr == 2 && println(io, " Maximum number of function evaluations (nfevalmax) reached.")
    end
    println(io, "")
    println(io, " Final objective function value = ", R.f)
    print(io, " Sample of best point = $(typeof(R.x))[ ")
    for i = 1:min(length(R.x) - 1, 3)
        print(io, "$(R.x[i]), ")
    end
    if length(R.x) > 4
        print(io, "..., ")
    end
    println(io, R.x[length(R.x)], "]")
    println(io, " Projected gradient norm = ", R.gnorm)
    println(io, "")
    println(io, " Number of iterations = ", R.nit)
    print(io, " Number of function evaluations = ", R.nfeval)
end
