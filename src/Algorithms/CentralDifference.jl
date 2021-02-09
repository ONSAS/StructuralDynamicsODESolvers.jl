"""
    CentralDifference{N} <: AbstractIntegrationAlgorithm

Central difference scheme with given step-size.

### Fields

- `Δt` -- step-size
"""
struct CentralDifference{N} <: AbstractIntegrationAlgorithm
    Δt::N
end

CentralDifference(; Δt::N) where N = CentralDifference(Δt)

function init(alg::CentralDifference, M, C, K, U₀, U₀′, U₀′′)
    Δt = alg.Δt

    # compute integration constants
    a₀ = 1/Δt^2
    a₁ = 1/(2*Δt)
    a₂ = 2*a₀
    a₃ = 1/a₂

    U⁻ = U₀ - Δt * U₀′ + a₃ * U₀′′
    M̂ = a₀*M + a₁*C

    return a₀, a₁, a₂, U⁻, M̂
end


"""
    solve(ivp::InitialValueProblem, alg, args..; kwargs...)

Solve an initial-value problem.
"""
function solve end

function solve(ivp::InitialValueProblem{<:SecondOrderAffineContinuousSystem{N}, XT},
               alg::CentralDifference{N},
               NSTEPS::Int) where {N, VT, XT<:Tuple{VT, VT}}

    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)
    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)
    a₀, a₁, a₂, U⁻, M̂ = init(alg, M, C, K, U₀, U₀′, U₀′′)
    M̂⁻¹ = inv(M̂)

    U = Vector{VT}(undef, IMAX+1)
    U[1] = U⁻
    U[2] = U₀

    @inbounds for i in 2:IMAX
        R̂ᵢ = R[i] - (K - a₂ * M) * U[i] - (a₀*M - a₁*C) * U[i-1]
        U[i+1] = M̂⁻¹ * R̂ᵢ
    end
    return IntegrationSolution(alg, view(U, 2:IMAX+1))
end
