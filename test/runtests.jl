usign Pkg
Pkg.add("Test")

using SPGBox
using Test

func(x) = x[1]^2 + x[2]^4

function grad!(x,g) 
  g[1] = 2*x
  g[2] = 4*x^3
end

x = [ 10. , 18. ]
R = spgbox!(x,func,grad!)
@test R.f ≈ zeros(2)

x = [ 10. , 18. ]
R = spgbox!(x,func,grad!,l=[2,-Inf])
@test R.f ≈ 4. 
@test R.x ≈ [2.,0.]

x = [ 10. , 18. ]
R = spgbox!(x,func,grad!,l=[-Inf,2])
@test R.f ≈ 16.
@test R.x ≈ [0.,2.]

x = [ 10. , 18. ]
R = spgbox!(x,func,grad!,u=[+Inf,-2])
@test R.f ≈ 4. 
@test R.x ≈ [-2.,0.]

x = [ 10. , 18. ]
R = spgbox!(x,func,grad!,u=[+Inf,-2])
@test R.f ≈ 16. 
@test R.x ≈ [0.,-2.]







