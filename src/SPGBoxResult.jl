"""

`SPGBoxResult`: Data structure that contains the results of the `spgbox` 
optimization method. 

  `x`: Vector containing final point. 
 
  `f`: Final objective function value.

  `gnorm`: Norm of the projected gradient on the bounds.

  `nit`: Number of iterations performed. 

  `nfeval`: Number of function evaluations.

  `ierr`: Status of execution: 
         0. Converged successfully to solution. 
         1. Maximum number of iterations achieved.
         2. Maximum number of function evaluations achieved.

"""
struct SPGBoxResult{T,F,G}
    x::T
    f::F
    gnorm::G
    nit::Int64
    nfeval::Int64
    ierr::Int64
end

function Base.show(io::IO, R::SPGBoxResult)
    println("")
    println(" SPGBOX RESULT: ")
    println("")
    R.ierr == 0 && println(" Convergence achieved. ")
    R.ierr == 1 && println(" Maximum number of iterations (nitmax) reached.")
    R.ierr == 2 && println(" Maximum number of function evaluations (nfevalmax) reached.")
    println("")
    println(" Final objective function value = ", R.f)
    print(" Sample of best point = $(typeof(R.x))[ ")
    for i = 1:min(length(R.x) - 1, 3)
        print("$(R.x[i]), ")
    end
    if length(R.x) > 4
        print("..., ")
    end
    println(R.x[length(R.x)], "]")
    println(" Projected gradient norm = ", R.gnorm)
    println("")
    println(" Number of iterations = ", R.nit)
    println(" Number of function evaluations = ", R.nfeval)
end
