using Test
using SPGBox
using Unitful
using ReverseDiff

@testset "simple quadratic" begin

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
    R = spgbox!(f, g!, x, lower=[2, -Inf])
    @test R.f ≈ 4.0
    @test R.x ≈ [2.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower=[-Inf, 2.0])
    @test R.f ≈ 1.0
    @test R.x ≈ [0.0, 2.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, upper=[+Inf, -2])
    @test R.f ≈ 9.0
    @test R.x ≈ [0.0, -2.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, upper=[-2, +Inf])
    @test R.f ≈ 4.0
    @test R.x ≈ [-2.0, 1.0]

    # just testing the interface
    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower=[-Inf, 2.0])
    @test R.f ≈ 1.0
    @test R.x ≈ [0.0, 2.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, upper=[-2, +Inf])
    @test R.f ≈ 4.0
    @test R.x ≈ [-2.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(f, g!, x, lower=[-Inf, 2.0], upper=[-2.0, +Inf])
    @test R.f ≈ 5.0
    @test R.x ≈ [-2.0, 2.0]

    # Test the mutating call  
    x = [10.0, 18.0]
    R = spgbox(f, g!, x, lower=[-Inf, -2.0], upper=[+Inf, 2])
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

    # Matrix input support
    x = [10.0 18.0]
    R = spgbox!(f, g!, x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0 1.0]

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
    x = [10.0, 18.0]u"nm"
    R = spgbox!(f_units, g_units!, x, lower=[-Inf, 2.0]u"nm")
    @test R.f ≈ 1.0u"nm^2"
    @test R.x ≈ [0.0, 2.0]u"nm"
    R = spgbox!(f_units, g_units!, x, lower=[-Inf, 2.0]u"nm")
    @test R.f ≈ 1.0u"nm^2"
    @test R.x ≈ [0.0, 2.0]u"nm"
    R = spgbox!(f_units, g_units!, x, upper=[-5.0, +Inf]u"nm")
    @test R.f ≈ 25.0u"nm^2"
    @test R.x ≈ [-5.0, 1.0]u"nm"
    R = spgbox!(f_units, g_units!, x, lower=[-Inf, 2.0]u"nm", upper=[-5.0, +Inf]u"nm")
    @test R.f ≈ 26.0u"nm^2"
    @test R.x ≈ [-5.0, 2.0]u"nm"

    #
    # With a single function to compute function and gradient
    #
    function fg!(g, x)
        g[1] = 2 * x[1]
        g[2] = 2 * (x[2] - 1)
        fx = x[1]^2 + (x[2] - 1)^2
        return fx
    end

    x = [10.0, 18.0]
    R = spgbox!(fg!, x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(fg!, x, lower=[-Inf, 2.0], upper=[-2.0, +Inf])
    @test R.f ≈ 5.0
    @test R.x ≈ [-2.0, 2.0]

    x = [10.0, 18.0]
    R = spgbox(fg!, x)
    @test R.f ≈ 0.0
    @test R.x ≈ [0.0, 1.0]

    x = [10.0, 18.0]
    R = spgbox!(fg!, x, lower=[-Inf, 2.0], upper=[-2.0, +Inf])
    @test R.f ≈ 5.0
    @test R.x ≈ [-2.0, 2.0]

    #
    # An example that verify that the line search was fixed
    #
    x0 = [i / 9 for i = 1:8]
    fH368(x) = -1.0 * sum(x .^ 2) * sum(x .^ 4) + sum(x .^ 3)^2
    function gH368!(g, x)
        ReverseDiff.gradient!(g, fH368, x)
    end
    lower = zeros(8)
    upper = ones(8)
    R = spgbox!(fH368, gH368!, x0, lower=lower, upper=upper)
    g = similar(x0)
    gH368!(g, R.x)
    @test SPGBox.pr_gradnorm(g, R.x, lower, upper) <= 1.0e-5

    # Test return from callback
    x = [10.0, 18.0]
    R = spgbox!(f, g!, x; callback = (R) -> R.nit <= 1 ? false : true)
    @test R.f ≈ 388.99999931805274
    @test R.nit == 2
    @test R.x ≈ [9.999999991234612, 17.99999998509884]

end
