"""
$(TYPEDEF)

Backward Euler's integration scheme with given step-size.

### Fields

- `Δt` -- step-size

"""
struct BackwardEuler{N} <: AbstractSolver
    Δt::N
end

BackwardEuler(; Δt::N) where N = BackwardEuler(Δt)

function _init(alg::BackwardEuler, M, C, K)
    Δt = alg.Δt
    K̂ = K * Δt + C
    return K̂
end

function _solve(alg::BackwardEuler{N},
                ivp::InitialValueProblem{<:SecondOrderAffineContinuousSystem{N}, XT},
                NSTEPS::Int) where {N, VT, XT<:Tuple{VT, VT}}

    sys = system(ivp)
    (U₀, _) = initial_state(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    iszero(M) || throw(ArgumentError("this method assumes that the system is of first order"))

    K̂ = _init(alg, M, C, K)
    K̂⁻¹ = inv(K̂)

    # initialize
    U = Vector{VT}(undef, IMAX)
    U[1] = U₀

    @inbounds for i in 1:NSTEPS
        R̂ᵢ₊₁ = R[i+1] * Δt + C*U[i]
        # solve
        U[i+1] = K̂⁻¹ * R̂ᵢ₊₁
    end

    return _build_solution(alg, U, NSTEPS)
end

function _build_solution(alg::BackwardEuler{N}, U, NSTEPS) where {N}
    t = range(zero(N), step=alg.Δt, length=NSTEPS)
    return Solution(alg, U, t)
end
