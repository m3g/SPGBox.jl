using TestItemRunner
@run_package_tests

@testitem "Aqua.test_all" begin
    import Aqua
    Aqua.test_all(SPGBox)
end

@testitem "simple quadratic" begin
    using Unitful
    using ReverseDiff

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
    R = spgbox!(f, g!, x; callback=(R) -> R.nit <= 1 ? false : true)
    @test R.f ≈ 388.99999931805274
    @test R.nit == 2
    @test R.x ≈ [9.999999991234612, 17.99999998509884]

    # Test definition of m by the user
    f2(x) = x[1]^2 + x[2]^2 + 1
    g2!(g, x) = g .= 2 .* x
    x0 = [5.0, 3.0]
    R = spgbox!(f2, g2!, x0; m = 5)
    @test R.f ≈ 1.0
end

@testitem "VAux" begin
    using SPGBox
    using BenchmarkTools
    f(x) = x[1]^2 + x[2]^2 + 1
    g!(g, x) = g .= 2 .* x
    x0 = [5.0, 3.0]
    vaux = SPGBox.VAux(x0, 0.0)
    # Test error when m and VAux are not compatible
    R = spgbox!(f, g!, x0; m = 10, vaux=vaux)
    @test R.f ≈ 1.0
    @test_throws DimensionMismatch spgbox!(f, g!, x0; m = 5, vaux=vaux)
    # Test allocations
    a = @ballocated spgbox!($f, $g!, $x0; vaux=$vaux) setup = (samples = 1; evals = 1)
    @test a == 0
    vaux = SPGBox.VAux(x0, 0.0; m=10)
    a = @ballocated spgbox!($f, $g!, $x0; vaux=$vaux) setup = (samples = 1; evals = 1)
    @test a == 0
end

@testitem "Optimization.jl interface" begin
    using SPGBox
    using Optimization
    using ForwardDiff

    f(x, p) = (x[1] - p[1])^2 + (x[2] - 1)^2
    optf = OptimizationFunction(f, AutoForwardDiff())

    # Unbounded
    prob = OptimizationProblem(optf, [10.0, 18.0], [3.0])
    sol = solve(prob, SPGBoxOpt())
    @test sol.u ≈ [3.0, 1.0]
    @test sol.objective ≈ 0.0 atol = 1e-8
    @test sol.retcode == ReturnCode.Success
    @test prob.u0 == [10.0, 18.0]

    # Bounded, and options passed to spgbox
    prob = OptimizationProblem(optf, [10.0, 18.0], [0.0]; lb=[2.0, -Inf], ub=[Inf, Inf])
    sol = solve(prob, SPGBoxOpt(eps=1e-8))
    @test sol.u ≈ [2.0, 1.0]
    @test sol.objective ≈ 4.0

    # Maximization
    optf_max = OptimizationFunction((x, p) -> -f(x, p), AutoForwardDiff())
    prob = OptimizationProblem(optf_max, [10.0, 18.0], [3.0]; sense=Optimization.MaxSense)
    @test solve(prob, SPGBoxOpt()).u ≈ [3.0, 1.0]
    # the callback receives the user's (maximized) objective, not its negation
    objs = Float64[]
    solve(prob, SPGBoxOpt(); callback=(state, loss) -> (push!(objs, loss); false))
    @test all(<=(0.0), objs)
    @test objs[1] > -f([10.0, 18.0], [3.0])

    # Limits, callback and init/solve! interface
    prob = OptimizationProblem(optf, [10.0, 18.0], [3.0])
    sol = solve(prob, SPGBoxOpt(); maxiters=1)
    @test sol.retcode == ReturnCode.MaxIters
    @test sol.stats.iterations == 1
    @test sol.original.nit == 2 # spgbox! counts the initial point as iteration 1
    iters = Int[]
    sol = solve(prob, SPGBoxOpt(); callback=(state, loss) -> (push!(iters, state.iter); true))
    @test sol.retcode == ReturnCode.Terminated
    @test iters == [1]
    @test sol.stats.iterations == 1
    sol = solve!(init(prob, SPGBoxOpt()))
    @test sol.retcode == ReturnCode.Success

    # Time limit and function-evaluation limit
    sol = solve(prob, SPGBoxOpt(); maxtime=0.0)
    @test sol.retcode == ReturnCode.MaxTime
    sol = solve(prob, SPGBoxOpt(nfevalmax=1))
    @test sol.original.ierr == 2

    @test_throws ArgumentError SPGBoxOpt(nitmax=10)
    @test sprint(show, SPGBoxOpt(eps=1e-8)) == "SPGBoxOpt(eps = 1.0e-8)"
end
