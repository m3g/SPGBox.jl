#
# Algorithm of:
#
# NONMONOTONE SPECTRAL PROJECTED GRADIENT METHODS ON CONVEX SETS
# ERNESTO G. BIRGIN, JOSÉ MARIO MARTÍNEZ, AND MARCOS RAYDAN
# SIAM J. O. PTIM. Vol. 10, No. 4, pp. 1196-1211
#
# Originally implemented by J. M. Martínez (IMECC - UNICAMP)
# Translated to Julia by L. Martínez (IQ-UNICAMP)
#
"""

`spgbox!(x :: Vector{Real}, func, grad; l :: Vector{Real}, u :: Vector{Real})`

Minimizes function `func` starting from initial point `x`, given
the function to compute the gradient, `grad`. 

Optional lower and upper box bounds can be provided using optional arguments `l` and `u`.

Returns a structure of type `SPGBoxResult`, containing the best solution found
in `x` and the final objective function in `f`.

# Examples
```jldocstest
julia> func(x) = x[1]^2 + x[2]^2

julia> function grad!(x,g)
         g[1] = 2*x[1]
         g[2] = 2*x[2]
       end
```
## Without bounds

```jldocstest
julia> x = rand(2)

julia> spgbox!(x,func,grad!)

 SPGBOX RESULT: 

 Convergence achieved. 

 Final objective function value = 0.0
 Best solution found = [ 0.0, 0.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 3
```

## With bounds

```jldocstest
julia> x = 2 .+ rand(2)

julia> spgbox!(x,func,grad!,l=[2.,-Inf])

 SPGBOX RESULT:

 Convergence achieved.

 Final objective function value = 4.0
 Best solution found = [ 2.0, 0.0]
 Projected gradient norm = 0.0

 Number of iterations = 2
 Number of function evaluations = 3
```
"""
function spgbox!(x :: AbstractVector{Float64}, func, grad;
                 l = [ -Inf for i in 1:length(x) ],
                 u = [ +Inf for i in 1:length(x) ],
                 g = similar(x),
                 xn = similar(x),
                 gn = similar(x),
                 nitmax :: Int64 = 100,
                 nfevalmax :: Int64 = 1000,
                 eps :: Float64 = 1.e-5,
                 iprint :: Int64 = 0,
                 m :: Int64 = 10,
                 fprev = zeros(m),
                 project_x0 = true)

  # Number of variables
  n = length(x)

  # Project initial point into bounds if set
  if project_x0
    for i in 1:length(x)
      x[i] = max(x[i],l[i])
      x[i] = min(x[i],u[i])
    end
  end

  # Objective function at initial point
  f = func(x)
  nfeval = 1
  grad!(x,g)

  tspg = 1.
  for i in 1:m
    fprev[i] = f
  end

  # Iteration counter
  nit = 0
  while nit < nitmax

    if iprint > 0
      println("----------------------------------------------------------- ")
      println(" Iteration: ", nit )
      println(" x = ", x[1], " ... ", x[n] )
      println(" Objective function value = ", f)
    end
    
    # Compute gradient norm
    gnorm = 0.
    for i in 1:n
      z = max(l[i], min(u[i], x[i]-g[i])) - x[i]
      gnorm = max(gnorm, abs(z))
    end
    if iprint > 0
      println(" ")
      println(" Norm of the projected gradient = ", gnorm)
      println(" Number of function evaluations = ", nfeval)
    end

    # Stopping criteria
    if gnorm <= eps
      ierr= 0
      return SPGBoxResult(x,f,gnorm,nit,nfeval,ierr)
    end
    if nfeval >= nfevalmax
      ierr= 2
      return SPGBoxResult(x,f,gnorm,nit,nfeval,ierr)
    end

    t = tspg
    fref = maximum(fprev)

    if iprint > 2
      println(" fref = ", fref)
      println(" fprev = ", fprev)
      println(" t = ", t)
    end

    fn = 1.1*fref
    while( fn > fref )
      for i in 1:n
        xn[i] = x[i] - t*g[i]
        xn[i] = max(l[i],min(xn[i],u[i]))
      end 
      if iprint > 2
        println(" xn = ", xn[1], xn[2] )
      end
      nfeval = nfeval + 1
      fn = func(xn)
      if iprint > 2
        println(" f[xn] = ", fn, " fref = ", fref )
      end
      # Maximum number of function evaluations achieved
      if nfeval > nfevalmax 
        ierr = 2
        return SPGBoxResult(x,f,gnorm,nit,nfeval,ierr)
      end
      # Reduce region
      t = t/2
    end

    grad!(xn,gn)
    num = 0.
    den = 0.
    for i in 1:n
      num = num + (xn[i]-x[i])^2
      den = den + (xn[i]-x[i])*(gn[i]-g[i])
    end
    if den <= 0.
      tspg = 100.
    else
      tspg =  min(1.e3,num/den)
    end
    f = fn
    for i in 1:n
      x[i] = xn[i]
      g[i] = gn[i]
    end
    for i in 1:m-1
      fprev[i] = fprev[i+1]
    end
    fprev[m] = f
    nit = nit + 1
  end

  # Maximum number of iterations achieved
  ierr = 1
  return SPGBoxResult(x,f,gnorm,nit,nfeval,ierr)

end

