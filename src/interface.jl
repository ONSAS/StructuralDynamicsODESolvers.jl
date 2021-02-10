"""
    AbstractSolver

Abstract supertype for all direct integration methods.
"""
abstract type AbstractSolver end

"""
    step_size(alg::AbstractSolver)

Return the step size of the given algorithm.

### Input

- `alg` -- algorithm

### Output

The step size of the algorithm, or `nothing` if the step-size is not fixed.
"""
step_size(alg::AbstractSolver) = nothing

"""
    AbstractSolution

Abstract supertype that holds the solution of a numerical integration.
"""
abstract type AbstractSolution end

"""
    Solution{T<:AbstractSolver, UT, VT, AT} <: AbstractSolution

### Fields

- `alg` -- algorithm used in the integration
- `U`   -- displacements
- `U′`  -- velocities
- `U′′` -- accelerations
- `t`   -- vector of time values
"""
struct Solution{T<:AbstractSolver, UT, VT, AT, ST} <: AbstractSolution
    alg::T
    U::UT
    U′::VT
    U′′::AT
    t::ST
end

# constructor with missing fields
Solution(alg, U, t) = Solution(alg, U, nothing, nothing, t)

"""
    displacements(sol::Solution)

Return the vector of displacements of the given solution.
"""
displacements(sol::Solution) = sol.U

"""
    velocities(sol::Solution)

Return the vector of velocities of the given solution.
"""
velocities(sol::Solution) = sol.U′

"""
    accelerations(sol::Solution)

Return the vector of accelerations of the given solution.
"""
accelerations(sol::Solution) = sol.U′′

"""
    times(sol::Solution)

Return the vector of times of the given solution.
"""
times(sol::Solution) = sol.t

# ===============
# Problem types
# ===============

abstract type AbstractProblem end

struct StructuralDynamicsProblem{T, PT} <: AbstractProblem
    alg::T
    ivp::PT
    NSTEPS::Int
end

"""
    solve(ivp::InitialValueProblem, alg::AbstractSolver, args..; kwargs...)

Solve an initial-value problem.

### Input

- `ivp` -- initial-value problem
- `alg` -- algorithm

### Output

A solution structure (`Solution`) that holds the result and the algorithm used
to obtain it.
"""
function solve(ivp::InitialValueProblem, alg::AbstractSolver, args...; NSTEPS, kwargs...)
    return solve!(init(ivp, alg, args...; NSTEPS=NSTEPS, kwargs...))
end

# internal defs
function solve!(prob::StructuralDynamicsProblem)
    return _solve(prob.alg, prob.ivp, prob.NSTEPS)
end

function init(ivp::InitialValueProblem{<:SecondOrderAffineContinuousSystem{N}, XT},
              alg::AbstractSolver;
              NSTEPS) where {N, VT, XT<:Tuple{VT, VT}}

    return StructuralDynamicsProblem(alg, ivp, NSTEPS)
end

# lazily extend the vector to the required number of steps
function _init_input(R::AbstractVector{N}, IMAX) where {N<:Number}
    return Fill(R, IMAX)
end

# no-op, checking that the number of forcing terms has the correct length
function _init_input(R::AbstractVector{VT}, IMAX) where {N, VT<:AbstractVector{N}}
    @assert length(R) == IMAX "expected the forcing term to be an array of length $IMAX, got $(length(R))"
    return R
end

# unwrap a second order system into its each component
function _unwrap(sys, IMAX)
    M = mass_matrix(sys)
    C = viscosity_matrix(sys)
    K = stiffness_matrix(sys)
    R = affine_term(sys)
    R = _init_input(R, IMAX)
    return M, C, K, R
end
