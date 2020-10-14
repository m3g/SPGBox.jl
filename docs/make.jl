import Pkg
Pkg.add("Documenter")
using Documenter
using SPGBox
push!(LOAD_PATH,"../src/")
makedocs(
    modules=[SPGBox],
    sitename="SPGBox.jl",
    pages = [
        "Home" => "index.md",
        "Installation" => "installation.md",
        "Basic usage" => "usage.md",
        "Options" => "options.md",
        "Reference" => "reference.md",
    ]
)
deploydocs(
    repo = "github.com/m3g/SPGBox.git",
    target = "build",
    branch = "gh-pages",
    versions = ["stable" => "v^", "v#.#" ],
)
