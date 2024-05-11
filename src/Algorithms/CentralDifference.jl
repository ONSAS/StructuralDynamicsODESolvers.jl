"""
$(TYPEDEF)

Central difference scheme with given step-size.

### Fields

- `Δt` -- step-size
"""
struct CentralDifference{N} <: AbstractSolver
    Δt::N
end

CentralDifference(; Δt::N) where {N} = CentralDifference(Δt)

function _init(alg::CentralDifference)
    Δt = alg.Δt

    # compute integration constants
    a₀ = 1 / Δt^2
    a₁ = 1 / (2 * Δt)
    a₂ = 2 * a₀
    a₃ = 1 / a₂

    return a₀, a₁, a₂, a₃
end

function _build_solution(alg::CentralDifference{N}, U, NSTEPS) where {N}
    t = range(zero(N); step=alg.Δt, length=(NSTEPS + 1))
    return Solution(alg, U, nothing, nothing, t)
end

function _solve(alg::CentralDifference{N},
                ivp::InitialValueProblem{ST,XT},
                NSTEPS::Int;
                kwargs...,) where {N,VT,ST,XT<:Tuple{VT,VT}}
    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)
    Δt = alg.Δt
    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)

    # integration constants
    a₀, a₁, a₂, a₃ = _init(alg)
    U⁻ = U₀ - Δt * U₀′ + a₃ * U₀′′
    M̂ = a₀ * M + a₁ * C
    M̂⁻¹ = factorize(M̂)
    M̂₋ = a₀ * M - a₁ * C

    U = Vector{VT}(undef, IMAX + 1)
    U[1] = U⁻
    U[2] = U₀

    @inbounds for i in 2:IMAX
        R̂ᵢ = R[i] - (K - a₂ * M) * U[i] - M̂₋ * U[i - 1]
        U[i + 1] = M̂⁻¹ \ R̂ᵢ
    end

    return _build_solution(alg, view(U, 2:(IMAX + 1)), NSTEPS)
end

function _solve(alg::CentralDifference{N},
                ivp::InitialValueProblem{ST,XT},
                NSTEPS::Int;
                t0=zero(N),
                kwargs...,) where {N,VT,ST<:SecondOrderContinuousSystem,XT<:Tuple{VT,VT}}
    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)
    Δt = alg.Δt
    IMAX = NSTEPS + 1
    @unpack M, C, fi, fe = sys

    # initial acceleration, obtained from the equilibrium
    # condition Mu'' + Cu' + fint(u) = fext(t), at time t = t0
    U₀′′ = M \ (fe(t0) - C * U₀′ - fi(U₀))

    # integration constants
    a₀, a₁, a₂, a₃ = _init(alg)
    U⁻ = U₀ - Δt * U₀′ + a₃ * U₀′′
    M̂ = a₀ * M + a₁ * C
    M̂⁻¹ = factorize(M̂)
    M̂₋ = a₀ * M - a₁ * C

    U = Vector{VT}(undef, IMAX + 1)
    U[1] = U⁻
    U[2] = U₀

    t = t0
    @inbounds for i in 2:IMAX
        R̂ᵢ = fe(t) - fi(U[i]) + a₂ * M * U[i] - M̂₋ * U[i - 1]
        U[i + 1] = M̂⁻¹ \ R̂ᵢ
        t += Δt
    end

    return _build_solution(alg, view(U, 2:(IMAX + 1)), NSTEPS)
end
