using Test
using SPGBox
using Unitful

@testset "simple polynomial" begin

    f(x) = x[1]^2 + (x[2] - 1)^2

    function g!(g, x)
        g[1] = 2 * x[1]
        g[2] = 2 * (x[2] - 1)
    end

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower = [2, -Inf])
    @test R.f ≈ 4.0
    @test R.x ≈ [2.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower = [-Inf, 2.0])
    @test R.f ≈ 1.0
    @test R.x ≈ [0.0, 2.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, upper = [+Inf, -2])
    @test R.f ≈ 9.0
    @test R.x ≈ [0.0, -2.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, upper = [-2, +Inf])
    @test R.f ≈ 4.0
    @test R.x ≈ [-2.0, 1.0]

    # just testing the interface
    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower = [-Inf, 2.0])
    @test R.f ≈ 1.0
    @test R.x ≈ [0.0, 2.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, upper = [-2, +Inf])
    @test R.f ≈ 4.0
    @test R.x ≈ [-2.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower = [-Inf, -2.0], upper = [+Inf, 2])
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, [-Inf, -2.0], [+Inf, 2], x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    # Test the mutating call  
    x = [10.0, 18.0]
    R = spgbox(f, g!, [-Inf, -2.0], [+Inf, 2], x)
    @test x == [10.0, 18.0]
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox(f, g!, x, lower = [-Inf, -2.0], upper = [+Inf, 2])
    @test x == [10.0, 18.0]
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    #BigFloat support
    x = BigFloat[10.0, 18.0]
    R = spgbox!(f, g!, x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]
    @test eltype(R.f) == eltype(R.x) == BigFloat

    #Float32 support
    x = Float32[10.0, 18.0]
    R = spgbox!(f, g!, x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]
    @test eltype(R.f) == eltype(R.x) == Float32

    # Unitful support
    f_units(x) = x[1]^2 + (x[2] - oneunit(eltype(x)))^2
    function g_units!(g, x)
        g[1] = 2 * x[1]
        g[2] = 2 * (x[2] - oneunit(eltype(x)))
    end
    x = [10.0, 18.0]u"nm"
    R = spgbox!(f_units, g_units!, x)
    @test R.f ≈ 0.0u"nm^2"
    @test R.x ≈ [0.0, 1.0]u"nm"

end
