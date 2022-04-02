using DataFrames
using CSV
using Plots
using BenchmarkProfiles

function genpprof(args)
    if length(args) == 1
        metric = args[1]
    else
        metric = "FEVAL"
    end
    base = DataFrame(CSV.File("base_bench.txt"))
    new = DataFrame(CSV.File("bench_results.txt"))

    pp = performance_profile(
    PlotsBackend(), 
    Float64.(hcat(base[!, metric], new[!, metric])),
    ["Base", "New"],
    title = "$metric profile"
    )
    savefig("profile.png")
end

if abspath(PROGRAM_FILE) == @__FILE__
    genpprof(ARGS)
end