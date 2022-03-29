using SPGBox
using NLPModels
using CUTEst
using CSV
using DataFrames
using BenchmarkTools
using Printf

function cutest2spg(probname)
    prob = CUTEstModel(probname)

    f(x) = obj(prob, x)
    function g!(g, x)
        return g .= grad(prob, x)
    end
    x0 = copy(prob.meta.x0)
    lower = copy(prob.meta.lvar)
    upper = copy(prob.meta.uvar)
    return f, g!, x0, lower, upper, prob
end

function cutest_test(probname)
    f, g!, x0, lower, upper, prob = cutest2spg(probname)
    R = spgbox!(f, g!, x0; lower=lower, upper=upper, nitmax=50_000, nfevalmax=200_000)
    stats = @benchmark spgbox!(
        $f, $g!, $x0; lower=$lower, upper=$upper, nitmax=50_000, nfevalmax=200_000
    )
    finalize(prob)
    return R, stats
end

function run_tests()
    tests = DataFrame(CSV.File("test_set.csv"))
    niters, nfevals, ns = Int[], Int[], Int[] 
    times = Float64[]

    for test in tests.TESTS
        println("Solving $test.")
        R, stats = cutest_test(test)
        solved = R.ierr == 0 ? 1 : -1
        append!(ns, length(R.x))
        append!(niters, solved * R.nit)
        append!(nfevals, solved * R.nfeval)
        append!(times, solved * median(stats.times) / 1.0e9)
    end
    results = copy(tests)
    results[!, "n"] .= ns
    results[!, "NITER"] .= niters
    results[!, "FEVAL"] .= nfevals
    results[!, "TIME"] .= times
    return results
end

function main()
    results = run_tests()
    # Format the time duration in scientific notation
    results[!, "TIME"] = map(t -> @sprintf("%e", t), results[!, "TIME"])
    CSV.write("bench_results.txt", results)
    return println(results)
end

main()
