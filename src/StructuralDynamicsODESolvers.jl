module StructuralDynamicsODESolvers

# dependencies
include("init.jl")

# User API
include("interface.jl")

# Available algorithms
include("Algorithms/Bathe.jl")
include("Algorithms/CentralDifference.jl")
include("Algorithms/Houbolt.jl")
include("Algorithms/Newmark.jl")
include("Algorithms/BackwardEuler.jl")

# exported methods and types
include("exports.jl")

end # module
