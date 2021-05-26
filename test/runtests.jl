using Test
using SPGBox

@testset "simple polynomial" begin

  f(x) = x[1]^2 + (x[2]-1.)^2
  
  function g!(g,x) 
    g[1] = 2*x[1]
    g[2] = 2*(x[2]-1.)
  end
  
  x = [ 10. , 18. ]
  R = spgbox!(f,g!,x)
  @test R.f ≈ 0.
  @test R.x ≈ [0.,1.]
  
  x = [ 10. , 18. ]
  R = spgbox!(f,g!,x,lower=[2,-Inf])
  @test R.f ≈ 4. 
  @test R.x ≈ [2.,1.]
  
  x = [ 10. , 18. ]
  R = spgbox!(f,g!,x,lower=[-Inf,2.])
  @test R.f ≈ 1.
  @test R.x ≈ [0.,2.]
  
  x = [ 10. , 18. ]
  R = spgbox!(f,g!,x,upper=[+Inf,-2])
  @test R.f ≈ 9. 
  @test R.x ≈ [0.,-2.]
  
  x = [ 10. , 18. ]
  R = spgbox!(f,g!,x,upper=[-2,+Inf])
  @test R.f ≈ 4. 
  @test R.x ≈ [-2.,1.]

end







