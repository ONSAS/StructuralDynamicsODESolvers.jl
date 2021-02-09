"""
    AbstractIntegrationAlgorithm

Abstract supertype for all direct integration methods.
"""
abstract type AbstractIntegrationAlgorithm end

"""
    step_size(alg::AbstractIntegrationAlgorithm)

Return the step size of the given algorithm.

### Input

- `alg` -- direct integration algorithm

### Output

The step size of the algorithm.
"""
step_size(alg::AbstractIntegrationAlgorithm) = alg.Δt

"""
    AbstractSolution
"""
abstract type AbstractSolution end

"""
    IntegrationSolution{T<:AbstractIntegrationAlgorithm, UT, VT, AT} <: AbstractSolution

### Fields
"""
struct IntegrationSolution{T<:AbstractIntegrationAlgorithm, UT, VT, AT} <: AbstractSolution
    alg::T
    U::UT
    U′::VT
    U′′::AT
end

# constructor with missing fields
IntegrationSolution(alg, U) = IntegrationSolution(alg, U, nothing, nothing)

"""
    displacements(sol::IntegrationSolution)
"""
displacements(sol::IntegrationSolution) = sol.U

"""
    velocities(sol::IntegrationSolution)
"""
velocities(sol::IntegrationSolution) = sol.U′

"""
    accelerations(sol::IntegrationSolution)
"""
accelerations(sol::IntegrationSolution) = sol.U′′

# lazily extend the vector to the required number of steps
function _init_input(R::AbstractVector{N}, IMAX) where {N<:Number}
    return Fill(R, IMAX)
end

# no-op, checking that the number of forcing terms has the correct length
function _init_input(R::AbstractVector{VT}, IMAX) where {N, VT<:AbstractVector{N}}
    @assert length(R) == IMAX "expected the forcing term to be an array of length $IMAX, got $(length(R))"
    return R
end

function _unwrap(sys, IMAX)
    M = mass_matrix(sys)
    C = viscosity_matrix(sys)
    K = stiffness_matrix(sys)
    R = affine_term(sys)
    R = _init_input(R, IMAX)
    return M, C, K, R
end
