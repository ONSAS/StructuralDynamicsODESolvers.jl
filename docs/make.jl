using Documenter, StructuralDynamicsODESolvers

DocMeta.setdocmeta!(StructuralDynamicsODESolvers, :DocTestSetup,
                    :(using StructuralDynamicsODESolvers); recursive=true)

# generate Literate documentation
include("generate.jl")

makedocs(; format=Documenter.HTML(; prettyurls=haskey(ENV, "GITHUB_ACTIONS"),  # disable for local builds
                                  collapselevel=1),
         sitename="StructuralDynamicsODESolvers.jl",
         pages=["Home" => "index.md",
                "Algorithms" => Any["First-order problems" => "lib/first_order.md",
                                    "Second-order problems" => "lib/second_order.md"],
                "Examples" => Any["Example (Ch.9 Bathe)" => "models/example_9_1_Bathe.md",
                                  "Spring-mass" => "models/massDashpotSpring.md",
                                  "Dynamic Von Mises" => "models/dynamic_von_mises_truss.md"],
                "API Reference" => "lib/api.md",
                "References" => "references.md",
                "About" => "about.md"])

deploydocs(; repo="github.com/ONSAS/StructuralDynamicsODESolvers.jl.git",
           push_preview=true)
