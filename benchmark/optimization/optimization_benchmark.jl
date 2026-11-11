#
# Benchmark of SPGBox (through its Optimization.jl interface, `SPGBoxOpt`)
# against other gradient-based, box-constrained local solvers:
#
# - `OptimizationLBFGSB.LBFGSB()`   : the classical L-BFGS-B algorithm.
# - `Optim.Fminbox(Optim.LBFGS())`  : L-BFGS with box constraints enforced by
#                                     an interior/barrier approach.
# - `NLopt.LD_LBFGS()`              : NLopt's L-BFGS, with bounds supported
#                                     directly by all NLopt local algorithms.
#
# These are the solvers most comparable to SPGBox: all are deterministic,
# gradient-based, local methods that natively support box constraints. 
# (SPGBox itself is a first-order method -- a spectral projected gradient --
# so it is expected to need more function/gradient evaluations than the
# quasi-Newton solvers, in exchange for a much lower cost per iteration and
# no need to store or update an approximate Hessian.)
#
# Run with:
#
#   julia --project=. optimization_benchmark.jl
#
# This generates `results.csv` with the results of every (problem, dimension,
# solver) combination, and prints a summary table to the terminal.
#
using SPGBox
using Optimization
using OptimizationLBFGSB
using OptimizationOptimJL
using OptimizationNLopt
using BenchmarkTools
using DataFrames
using CSV
using Printf

#
# Test problems: classic scalable test functions, with box bounds that force
# part of the solution onto the boundary (except for the Powell singular
# function, for which the bounds are inactive at the solution, to test
# behavior on an ill-conditioned problem without boundary effects).
#
# Each problem returns: the initial point `x0`, and the `lower` and `upper`
# bound vectors (all feasible with respect to each other).
#

function rosenbrock(x, p)
    n = length(x)
    s = zero(eltype(x))
    for i in 1:2:n-1
        s += 100 * (x[i+1] - x[i]^2)^2 + (1 - x[i])^2
    end
    return s
end
function rosenbrock_grad!(g, x, p)
    n = length(x)
    fill!(g, 0.0)
    for i in 1:2:n-1
        g[i] += -400 * x[i] * (x[i+1] - x[i]^2) - 2 * (1 - x[i])
        g[i+1] += 200 * (x[i+1] - x[i]^2)
    end
    return g
end
function rosenbrock_setup(n)
    iseven(n) || error("rosenbrock requires an even n")
    x0 = [isodd(i) ? -1.2 : 0.5 for i in 1:n] # unconstrained optimum is 1.0 for all i
    lower = fill(-5.0, n)
    upper = fill(0.8, n) # forces the even-indexed components onto the boundary
    return x0, lower, upper
end

function powell_singular(x, p)
    n = length(x)
    s = zero(eltype(x))
    for j in 1:4:n-3
        x1, x2, x3, x4 = x[j], x[j+1], x[j+2], x[j+3]
        f1 = x1 + 10x2
        f2 = sqrt(5) * (x3 - x4)
        f3 = (x2 - 2x3)^2
        f4 = sqrt(10) * (x1 - x4)^2
        s += f1^2 + f2^2 + f3^2 + f4^2
    end
    return s
end
function powell_singular_grad!(g, x, p)
    n = length(x)
    for j in 1:4:n-3
        x1, x2, x3, x4 = x[j], x[j+1], x[j+2], x[j+3]
        g[j] = 2 * (x1 + 10x2) + 40 * (x1 - x4)^3
        g[j+1] = 20 * (x1 + 10x2) + 4 * (x2 - 2x3)^3
        g[j+2] = 10 * (x3 - x4) - 8 * (x2 - 2x3)^3
        g[j+3] = -10 * (x3 - x4) - 40 * (x1 - x4)^3
    end
    return g
end
function powell_singular_setup(n)
    n % 4 == 0 || error("powell_singular requires n to be a multiple of 4")
    # half of the standard starting point (3,-1,0,1,...), scaled to be
    # feasible with the (inactive at the optimum, x=0) bounds below
    x0 = [[1.5, -0.5, 0.0, 0.5][mod1(i, 4)] for i in 1:n]
    lower = fill(-2.0, n)
    upper = fill(2.0, n)
    return x0, lower, upper
end

function diagquad(x, p)
    s = zero(eltype(x))
    for i in eachindex(x)
        s += i * (x[i] - 1)^2
    end
    return s
end
function diagquad_grad!(g, x, p)
    for i in eachindex(x)
        g[i] = 2i * (x[i] - 1)
    end
    return g
end
function diagquad_setup(n)
    x0 = zeros(n) # unconstrained optimum is 1.0 for all i
    lower = fill(-5.0, n)
    upper = [i <= n ÷ 2 ? 0.5 : 5.0 for i in 1:n] # forces the first half onto the boundary
    return x0, lower, upper
end

const problems = [
    ("Rosenbrock", rosenbrock, rosenbrock_grad!, rosenbrock_setup),
    ("Powell singular", powell_singular, powell_singular_grad!, powell_singular_setup),
    ("Diagonal quadratic", diagquad, diagquad_grad!, diagquad_setup),
]

const solvers = [
    ("SPGBox", SPGBoxOpt()),
    ("LBFGSB", OptimizationLBFGSB.LBFGSB()),
    ("Fminbox(LBFGS)", Optim.Fminbox(Optim.LBFGS())),
    ("NLopt.LD_LBFGS", NLopt.LD_LBFGS()),
]

const dims = [20, 100, 1000]

#
# For each (problem, dimension, solver), solve the problem and record the
# final objective value, the number of function and gradient evaluations
# (counted independently of each solver's own reporting, since not all
# solvers report evaluation counts through `sol.stats`), and the median
# solve time.
#
function run_benchmarks(; maxiters=5000, samples=5)
    results = DataFrame(
        problem=String[], n=Int[], solver=String[],
        f=Float64[], retcode=String[], nf=Int[], ng=Int[], time=Float64[],
    )
    for (pname, f, g!, setup) in problems, n in dims
        x0, lower, upper = setup(n)
        for (sname, alg) in solvers
            nf = Ref(0)
            ng = Ref(0)
            f_counted(x, p) = (nf[] += 1; f(x, p))
            g_counted!(g, x, p) = (ng[] += 1; g!(g, x, p))
            optf = OptimizationFunction(f_counted; grad=g_counted!)
            prob = OptimizationProblem(optf, x0, nothing; lb=lower, ub=upper)
            local sol
            try
                sol = solve(prob, alg; maxiters=maxiters)
            catch e
                @warn "failed" problem = pname n solver = sname exception = e
                continue
            end
            t = @belapsed solve($prob, $alg; maxiters=$maxiters) evals = 1 samples = samples
            nf[] = 0
            ng[] = 0
            solve(prob, alg; maxiters=maxiters) # re-run once to recover accurate nf/ng
            push!(results, (
                problem=pname, n=n, solver=sname,
                f=sol.objective, retcode=string(sol.retcode),
                nf=nf[], ng=ng[], time=t,
            ))
            @printf(
                "%-20s n=%-5d %-16s f=%-12.6g retcode=%-10s nf=%-6d ng=%-6d t=%.4gs\n",
                pname, n, sname, sol.objective, string(sol.retcode), nf[], ng[], t,
            )
        end
    end
    return results
end

function benchmark()
    results = run_benchmarks()
    CSV.write("results.csv", results)
    return println(results)
end

if abspath(PROGRAM_FILE) == @__FILE__
    benchmark()
end
