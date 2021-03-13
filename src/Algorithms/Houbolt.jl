"""
$(TYPEDEF)

Houbolt's integration scheme with given step-size.

### Fields

- `Δt` -- step-size

### References

See [[HOU50]](@ref).
"""
struct Houbolt{N} <: AbstractSolver
    Δt::N
end

Houbolt(; Δt::N) where N = Houbolt(Δt)

function _init(alg::Houbolt, M, C, K)
    Δt = alg.Δt

    # compute integration constants
    a₀ = 2/Δt^2
    a₁ = 11/(6*Δt)
    a₂ = 5/Δt^2
    a₃ = 3/Δt
    a₄ = -2*a₀
    a₅ = -a₃/2
    a₆ = a₀/2
    a₇ = a₃/9

    K̂ = K + a₀ * M + a₁ * C

    return a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂
end

function _solve(alg::Houbolt{N},
                ivp::InitialValueProblem{<:SecondOrderAffineContinuousSystem{N}, XT},
                NSTEPS::Int) where {N, VT, XT<:Tuple{VT, VT}}

    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)
    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂ = _init(alg, M, C, K)
    K̂⁻¹ = inv(K̂)

    U = Vector{VT}(undef, IMAX)
    U[1] = U₀

    # special starting procedure to calculate U(Δt) and U(2Δt)
    V = _solve(CentralDifference(Δt=alg.Δt), ivp, 2) |> displacements
    U[2] = V[2]
    U[3] = V[3]

    @inbounds for i in 3:NSTEPS
        mᵢ = M * (a₂ * U[i] + a₄ * U[i-1] + a₆ * U[i-2])
        cᵢ = C * (a₃ * U[i] + a₅ * U[i-1] + a₇ * U[i-2])
        R̂ᵢ₊₁ = R[i+1] + mᵢ + cᵢ
        U[i+1] = K̂⁻¹ * R̂ᵢ₊₁
    end

    return _build_solution(alg, U, NSTEPS)
end

function _build_solution(alg::Houbolt{N}, U, NSTEPS) where {N}
    t = range(zero(N), step=alg.Δt, length=(NSTEPS+1))
    return Solution(alg, U, nothing, nothing, t)
end
