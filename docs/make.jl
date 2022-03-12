import Pkg
Pkg.add("Documenter")
using Documenter
using SPGBox
push!(LOAD_PATH, "../src/")
makedocs(
    modules = [SPGBox],
    sitename = "SPGBox.jl",
    pages = [
        "Installation" => "index.md",
        "User guide" => "usage.md",
        "Options" => "options.md",
        "Reference" => "reference.md",
    ],
)
deploydocs(
    repo = "github.com/m3g/SPGBox.jl.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#"],
)
