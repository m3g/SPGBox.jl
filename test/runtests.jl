import Pkg
Pkg.add("Test")
Pkg.add("SPGBox")
using Test
using SPGBox

@testset "simple polynomial" begin

  func(x) = x[1]^2 + (x[2]-1.)^2
  
  function grad!(x,g) 
    g[1] = 2*x[1]
    g[2] = 2*(x[2]-1.)
  end
  
  x = [ 10. , 18. ]
  R = spgbox!(x,func,grad!)
  @test R.f ≈ 0.
  @test R.x ≈ [0.,1.]
  
  x = [ 10. , 18. ]
  R = spgbox!(x,func,grad!,l=[2,-Inf])
  @test R.f ≈ 4. 
  @test R.x ≈ [2.,1.]
  
  x = [ 10. , 18. ]
  R = spgbox!(x,func,grad!,l=[-Inf,2.])
  @test R.f ≈ 1.
  @test R.x ≈ [0.,2.]
  
  x = [ 10. , 18. ]
  R = spgbox!(x,func,grad!,u=[+Inf,-2])
  @test R.f ≈ 9. 
  @test R.x ≈ [0.,-2.]
  
  x = [ 10. , 18. ]
  R = spgbox!(x,func,grad!,u=[-2,+Inf])
  @test R.f ≈ 4. 
  @test R.x ≈ [-2.,1.]

end







