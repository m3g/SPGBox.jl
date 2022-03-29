using DataFrames
using CSV
using Plots
using BenchmarkProfiles

base = DataFrame(CSV.File("base_bench.txt"))
new = DataFrame(CSV.File("bench_results.txt"))
pp = performance_profile(
    PlotsBackend(), 
    Float64.(hcat(base[!, "FEVAL"], new[!, "FEVAL"])),
    ["Base", "New"],
    title = "Simple comparison"
)
savefig("pp.png")
