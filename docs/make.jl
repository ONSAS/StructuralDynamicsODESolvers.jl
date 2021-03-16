ENV["GKSwstype"] = "100"  # set 'GR environment' to 'no output' (for Travis CI)
using Documenter, StructuralDynamicsODESolvers

DocMeta.setdocmeta!(StructuralDynamicsODESolvers, :DocTestSetup,
                   :(using StructuralDynamicsODESolvers); recursive=true)

# generate Literate documentation
include("generate.jl")

makedocs(
    format = Documenter.HTML(prettyurls = haskey(ENV, "GITHUB_ACTIONS"),  # disable for local builds
                             collapselevel = 1),
    sitename = "StructuralDynamicsODESolvers.jl",
    doctest = false,
    strict = false,
    pages = [
        "Home" => "index.md",
        "Algorithms" => Any["First-order problems" => "lib/first_order.md",
                            "Second-order problems" => "lib/second_order.md"],
        "Examples" => Any["Example (Ch.9 Bathe)" => "models/example_9_1_Bathe.md",
                          "Spring-mass" => "models/massDashpotSpring.md"],
        "API Reference" => "lib/api.md",
        "References" => "references.md",
        "About" => "about.md"
    ]
)

# Deploy built documentation from Travis.
deploydocs(
    repo = "github.com/ONSAS/StructuralDynamicsODESolvers.jl.git",
    push_preview = true,
)
