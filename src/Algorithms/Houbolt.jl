"""
    Houbolt{N} <: AbstractIntegrationAlgorithm

Houbolt's integration scheme with given step-size.

### Fields

- `Δt` -- step-size
"""
struct Houbolt{N} <: AbstractIntegrationAlgorithm
    Δt::N
end

Houbolt(; Δt::N) where N = Houbolt(Δt)

function init(alg::Houbolt, M, C, K)
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

function solve(ivp::InitialValueProblem{<:SecondOrderAffineContinuousSystem{N}, XT},
               alg::Houbolt{N},
               NSTEPS::Int) where {N, VT, XT<:Tuple{VT, VT}}

    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)
    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇, K̂ = init(alg, M, C, K)
    K̂⁻¹ = inv(K̂)

    U = Vector{VT}(undef, IMAX)
    U[1] = U₀

    # special starting procedure to calculate U(Δt) and U(2Δt)
    V = solve(ivp, CentralDifference(Δt=alg.Δt), 2) |> displacements
    U[2] = V[2]
    U[3] = V[3]

    @inbounds for i in 3:NSTEPS
        mᵢ = M * (a₂ * U[i] + a₄ * U[i-1] + a₆ * U[i-2])
        cᵢ = C * (a₃ * U[i] + a₅ * U[i-1] + a₇ * U[i-2])
        R̂ᵢ₊₁ = R[i+1] + mᵢ + cᵢ
        U[i+1] = K̂⁻¹ * R̂ᵢ₊₁
    end
    return IntegrationSolution(alg, U)
end
