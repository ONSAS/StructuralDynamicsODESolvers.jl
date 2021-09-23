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
    dim(sol::Solution)

Return the ambient dimension of the state space of the solution.
"""
dim(sol::Solution) = length(first(sol.U))

"""
    displacements(sol::Solution)

Return the vector of displacements of the given solution.
"""
displacements(sol::Solution) = sol.U

"""
    displacements(sol::Solution, i::Int)

Return the vector of displacements of the given solution along coordinate `i`.
"""
function displacements(sol::Solution, i::Int)
    1 ≤ i ≤ dim(sol) || throw(ArgumentError("expected the coordinate to be between 1 and $(dim(sol)), got $i"))
    U = displacements(sol)
    return [u[i] for u in U]
end

"""
    velocities(sol::Solution)

Return the vector of velocities of the given solution.
"""
velocities(sol::Solution) = sol.U′

"""
    velocities(sol::Solution, i::Int)

Return the vector of velocities of the given solution along coordinate `i`.
"""
function velocities(sol::Solution, i::Int)
    1 ≤ i ≤ dim(sol) || throw(ArgumentError("expected the coordinate to be between 1 and $(dim(sol)), got $i"))
    U′ = velocities(sol)
    return [u′[i] for u′ in U′]
end

"""
    accelerations(sol::Solution)

Return the vector of accelerations of the given solution.
"""
accelerations(sol::Solution) = sol.U′′

"""
    accelerations(sol::Solution, i::Int)

Return the vector of accelerations of the given solution along coordinate `i`.
"""
function accelerations(sol::Solution, i::Int)
    1 ≤ i ≤ dim(sol) || throw(ArgumentError("expected the coordinate to be between 1 and $(dim(sol)), got $i"))
    U′′ = accelerations(sol)
    return [u′′[i] for u′′ in U′′]
end

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
function solve(ivp::IVP{<:AbstractContinuousSystem}, alg::AbstractSolver, args...; NSTEPS, kwargs...)
    return solve!(init(ivp, alg, args...; NSTEPS=NSTEPS, kwargs...))
end

# internal defs
function solve!(prob::StructuralDynamicsProblem)
    return _solve(prob.alg, prob.ivp, prob.NSTEPS)
end

const SOACS = SecondOrderConstrainedLinearControlContinuousSystem
const SOCLCCS  = SecondOrderConstrainedLinearControlContinuousSystem

function init(ivp::InitialValueProblem{ST, XT},
              alg::AbstractSolver;
              NSTEPS) where {N, VT,
                             ST, # FIXME restrict to SOACS and SOCLCCS
                             XT<:Tuple{VT, VT}}

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

# TODO dispatch on B (eg. IdentityMultiple)
function _init_input(R::AbstractVector{VT}, B::AbstractMatrix, IMAX) where {N, VT<:AbstractVector{N}}
    @assert length(R) == IMAX "expected the forcing term to be an array of length $IMAX, got $(length(R))"
    return [B*Ri for Ri in R]
end

# unwrap a second order system into its each component
function _unwrap(sys::SecondOrderAffineContinuousSystem, IMAX)
    M = mass_matrix(sys)
    C = viscosity_matrix(sys)
    K = stiffness_matrix(sys)
    R = affine_term(sys)
    R = _init_input(R, IMAX)
    return M, C, K, R
end

function _unwrap(sys::SecondOrderConstrainedLinearControlContinuousSystem, IMAX)
    M = mass_matrix(sys)
    C = viscosity_matrix(sys)
    K = stiffness_matrix(sys)
    R = inputset(sys)
    B = input_matrix(sys)
    R = _init_input(R, B, IMAX)
    return M, C, K, R
end

# ================
# Plot recipes
# ================

function _check_vars(vars)
    if vars == nothing
        throw(ArgumentError("default ploting variables not implemented yet; you need " *
              "to pass the `vars=(...)` option, e.g. `vars=(0, 1)` to plot variable with " *
              "index 1 vs. time, or `vars=(1, 2)` to plot variable with index 2 vs. variable with index 1`"))
    end
    D = length(vars)
    @assert (D == 1) || (D == 2) "can only plot in one or two dimensions, " *
                                 "but received $D variable indices where `vars = ` $vars"
end

# plot displacements of the solution for the given vars tuple, eg. vars=(0, 1) for x1(t) vs t
@recipe function plot_solution(sol::Solution; vars=nothing, func=displacements)

   seriestype -->  :path # :scatter
   markershape --> :circle

   _check_vars(vars)

   if vars[1] == 0 && vars[2] != 0
       x = times(sol)
       y = func(sol, vars[2])
       x, y
    else
       x = func(sol, vars[1])
       y = func(sol, vars[2])
    end
    return x, y
end
