"""
$(TYPEDEF)

Bathe's integration scheme with sub-steps of equal size.

### Fields

- `Δt` -- step-size
- `α`  -- parameter α of the method
- `δ`  -- parameter δ of the method

### References

See [[BAT07]](@ref).
"""
struct Bathe{N} <: AbstractSolver
    Δt::N
end

Bathe(; Δt::N) where {N} = Bathe(Δt)

step_size(alg::Bathe) = alg.Δt

function _init(alg::Bathe, M, C, K)
    Δt = alg.Δt

    # compute integration constants
    ξ = 1/Δt
    ξ² = ξ*ξ
    a₀ = 16*ξ²
    a₁ = 4*ξ
    a₂ = 9*ξ²
    a₃ = 3*ξ
    a₄ = 2*a₁
    a₅ = 12*ξ²
    a₆ = -3*ξ²
    a₇ = -ξ

    K̂₁ = K + a₀ * M + a₁ * C
    K̂₂ = K + a₂ * M + a₃ * C
    return a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂₁, K̂₂
end

# a comment about notation: variables with superscript `+` correspond to
# evaluation at intermediate times, i.e. t + Δt/2
function _solve(alg::Bathe{N},
                ivp::InitialValueProblem{ST, XT},
                NSTEPS::Int) where {N, VT, ST, XT<:Tuple{VT, VT}}

    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)
    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂₁, K̂₂ = _init(alg, M, C, K)
    K̂₁⁻¹ = inv(K̂₁)
    K̂₂⁻¹ = inv(K̂₂)

    # initialize displacements, velocities and accelerations
    U = Vector{VT}(undef, IMAX)
    U′ = Vector{VT}(undef, IMAX)
    U′′ = Vector{VT}(undef, IMAX)
    U[1] = U₀
    U′[1] = U₀′
    U′′[1] = U₀′′

    @inbounds for i in 1:NSTEPS
        # ---------------------------------------------
        # First sub-step (Newmark trapezoidal rule)
        # ---------------------------------------------

        # calculate effective loads at time t + Δt/2
        mᵢ = M * (a₀ * U[i] + a₄ * U′[i] + U′′[i])
        cᵢ = C * (a₁ * U[i] + U′[i])
        Rᵢ⁺ = (R[i] + R[i+1])/2
        R̂ᵢ⁺ = Rᵢ⁺ + mᵢ + cᵢ

        # solve for displacements at time t + Δt/2
        Uᵢ⁺ = K̂₁⁻¹ * R̂ᵢ⁺

        # calculate velocities at time t + Δt/2
        U′ᵢ⁺ = a₁ * (Uᵢ⁺ - U[i]) - U′[i]

        # ---------------------------------------------
        # Second sub-step (3-point Euler backward method)
        # ---------------------------------------------

        # calculate effective loads at time t + Δt
        mᵢ = M * (a₅ * Uᵢ⁺ + a₆ * U[i] + a₁ * U′ᵢ⁺ + a₇ * U′[i])
        cᵢ = C * (a₁ * Uᵢ⁺ + a₇ * U[i])
        R̂ᵢ₊₁ = R[i+1] + mᵢ + cᵢ

        # solve for displacements at time t + Δt
        U[i+1] = K̂₂⁻¹ * R̂ᵢ₊₁

        # calculate velocities and accelerations at time t + Δt
        U′[i+1] = -a₇ * U[i] - a₁ * Uᵢ⁺ + a₃ * U[i+1]
        U′′[i+1] = -a₇ * U′[i] - a₁ * U′ᵢ⁺ + a₃ * U′[i+1]
    end

    return _build_solution(alg, U, U′, U′′, NSTEPS)
end

function _build_solution(alg::Bathe{N}, U, U′, U′′, NSTEPS) where {N}
    t = range(zero(N), step=alg.Δt, length=(NSTEPS+1))
    return Solution(alg, U, U′, U′′, t)
end
