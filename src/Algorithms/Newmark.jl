"""
$(TYPEDEF)

Newmark's integration scheme with given step-size and parameters `α` and `δ`.

### Fields

- `Δt` -- step-size
- `α`  -- parameter α of the method
- `δ`  -- parameter δ of the method

### References

See [[NEW59]](@ref).
"""
struct Newmark{N} <: AbstractSolver
    Δt::N
    α::N
    δ::N
end

function Newmark(; Δt::T, α::A, δ::D, check=true) where {T, A, D}
    N = promote_type(T, A, D)

    # consistency checks
    if check
        δ ≥ 0.5 || throw(ArgumentError("expected the parameter δ to be at least 1/2, but it is $δ"))
        μ = (1/2 + δ)^2 / 4
        α ≥ μ || throw(ArgumentError("expected the parameter α to be at least $μ, but it is $α"))
    end
    return Newmark(N(Δt), N(α), N(δ))
end

# notable special cases

"""
    Linear(Δt::N)

Linear integration scheme. Special case of Newmark with ``δ=1/2`` and ``α=1/6``.
"""
Linear(Δt::N) where N = Newmark(Δt=Δt, δ=1/2, α=1/6, check=false)
Linear(; Δt::N) where N = Linear(Δt)

"""
    Trapezoidal(Δt::N)

Trapezoidal integration scheme. Special case of Newmark with ``δ=1/2`` and ``α=1/4``.
"""
Trapezoidal(Δt::N) where N = Newmark(Δt=Δt, δ=1/2, α=1/4)
Trapezoidal(; Δt::N) where N = Trapezoidal(Δt)

function _init(alg::Newmark, M, C, K)
    Δt = alg.Δt
    α = alg.α
    δ = alg.δ

    # compute integration constants
    ξ = 1 / (α * Δt)
    a₀ = ξ/Δt
    a₁ = δ*ξ
    a₂ = ξ
    a₃ = 1/(2α) - 1
    a₄ = δ/α - 1
    a₅ = (Δt/2) * (a₄-1)
    a₆ = Δt*(1-δ)
    a₇ = δ*Δt

    K̂ = K + a₀ * M + a₁ * C

    return a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂
end

function _solve(alg::Newmark{N},
                ivp::InitialValueProblem{ST, XT},
                NSTEPS::Int) where {N, VT, ST, XT<:Tuple{VT, VT}}

    sys = system(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂ = _init(alg, M, C, K)
    K̂⁻¹ = factorize(K̂)

    # initialize displacements, velocities and accelerations
    (U₀, U₀′) = initial_state(ivp)
    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)

    U = Vector{VT}(undef, IMAX)
    U′ = Vector{VT}(undef, IMAX)
    U′′ = Vector{VT}(undef, IMAX)
    U[1] = U₀
    U′[1] = U₀′
    U′′[1] = U₀′′

    @inbounds for i in 1:NSTEPS
        # calculate effective loads
        mᵢ = M * (a₀ * U[i] + a₂ * U′[i] + a₃ * U′′[i])
        cᵢ = C * (a₁ * U[i] + a₄ * U′[i] + a₅ * U′′[i])
        R̂ᵢ₊₁ = R[i+1] + mᵢ + cᵢ

        # solve for displacements
        U[i+1] = K̂⁻¹ \ R̂ᵢ₊₁

        # calculate accelerations and velocities
        U′′[i+1] = a₀ * (U[i+1] - U[i]) - a₂ * U′[i] - a₃ * U′′[i]
        U′[i+1] = U′[i] + a₆ * U′′[i] + a₇ * U′′[i+1]
    end

    return _build_solution(alg, U, U′, U′′, NSTEPS)
end

function _build_solution(alg::Newmark{N}, U, U′, U′′, NSTEPS) where {N}
    t = range(zero(N), step=alg.Δt, length=(NSTEPS+1))
    return Solution(alg, U, U′, U′′, t)
end

function _solve_statistics(alg::Newmark{N},
                           ivp::InitialValueProblem{ST, XT},
                           NSTEPS::Int,
                           remake, idx) where {N, VT, ST, XT<:Tuple{VT, VT}}

    sys = system(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)
    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂ = _init(alg, M, C, K)
    K̂⁻¹ = factorize(K̂)

    # initialize displacements, velocities and accelerations
    U = Vector{VT}(undef, IMAX)
    U′ = Vector{VT}(undef, IMAX)
    U′′ = Vector{VT}(undef, IMAX)

    # output vectors
    U_min = fill(Inf, IMAX)
    U_max = fill(-Inf, IMAX)

    U′_min = fill(Inf, IMAX)
    U′_max = fill(-Inf, IMAX)

    for state in remake
        (U₀, U₀′) = state

        # initialization
        U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)
        U[1] = U₀
        U′[1] = U₀′
        U′′[1] = U₀′′

        @inbounds for i in 1:NSTEPS
            # calculate effective loads
            mᵢ = M * (a₀ * U[i] + a₂ * U′[i] + a₃ * U′′[i])
            cᵢ = C * (a₁ * U[i] + a₄ * U′[i] + a₅ * U′′[i])
            R̂ᵢ₊₁ = R[i+1] + mᵢ + cᵢ

            # solve for displacements
            U[i+1] = K̂⁻¹ \ R̂ᵢ₊₁

            # calculate accelerations and velocities
            U′′[i+1] = a₀ * (U[i+1] - U[i]) - a₂ * U′[i] - a₃ * U′′[i]
            U′[i+1] = U′[i] + a₆ * U′′[i] + a₇ * U′′[i+1]
        end

        # store observables
        @inbounds for i in 1:IMAX
            U_min[i] = min(U_min[i], U[i][idx])
            U_max[i] = max(U_max[i], U[i][idx])
            U′_min[i] = min(U′_min[i], U′[i][idx])
            U′_max[i] = max(U′_max[i], U′[i][idx])
        end
    end

    t = range(zero(N), step=alg.Δt, length=(NSTEPS+1))
    return SolutionExtrema(alg, U_min, U_max, U′_min, U′_max, idx, t)
end
