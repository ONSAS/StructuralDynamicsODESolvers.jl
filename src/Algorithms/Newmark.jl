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

function _build_solution(alg::Newmark{N}, U, U′, U′′, NSTEPS) where {N}
    t = range(zero(N), step=alg.Δt, length=(NSTEPS+1))
    return Solution(alg, U, U′, U′′, t)
end

function _init(alg::Newmark)
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

    return a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇
end

function _solve(alg::Newmark{N},
                ivp::InitialValueProblem{ST, XT},
                NSTEPS::Int) where {N, VT, ST, XT<:Tuple{VT, VT}}

    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)

    IMAX = NSTEPS + 1
    M, C, K, R = _unwrap(sys, IMAX)

    U₀′′ = M \ (R[1] - C * U₀′ - K * U₀)
    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇ = _init(alg)

    K̂ = K + a₀ * M + a₁ * C
    K̂⁻¹ = factorize(K̂)

    # initialize displacements, velocities and accelerations
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

# case with possibly non-constant fint(x) term
function _solve(alg::Newmark{N},
                ivp::InitialValueProblem{ST, XT},
                NSTEPS::Int,
                reltol=1e-6,
                maxiter=10,
                t0=zero(N)) where {N, VT, ST<:SecondOrderContinuousSystem, XT<:Tuple{VT, VT}}

    println("reltol = $reltol")
#function solve(problema::ProblemaNoLineal, alg::Newmark,
#               t0, tf, u0, v0;
#               reltol = 1e-6,  # tolerancia relativa en los desplazamientos
#               maxiter = 10)   # maximo numero de iteraciones

    sys = system(ivp)
    (U₀, U₀′) = initial_state(ivp)

    IMAX = NSTEPS + 1
    @unpack M, C, Fint, Fext = sys

    # initial acceleration, obtained from the equilibrium
    # condition Mu'' + Cu' + fint(u) = fext(t), at time t = t0
    U₀′′ = M \ (Fext(t0) - C * U₀′ - Fint(U₀))

    # integration constants
    a₀, a₁, a₂, a₃, a₄, a₅, a₆, a₇ = _init(alg)

    # initialize displacements, velocities and accelerations
    U = Vector{VT}(undef, IMAX)
    U′ = Vector{VT}(undef, IMAX)
    U′′ = Vector{VT}(undef, IMAX)
    U[1] = U₀
    U′[1] = U₀′
    U′′[1] = U₀′′

    # tangent matrix function
    KT(x) = ForwardDiff.jacobian(Fint, x)

    # times vector
    tvec = Vector{N}(undef, IMAX)
    tvec[1] = t0

    # iteration counter
    i = 1

    while i <= NSTEPS
        j = 1
        Ferr = +Inf
        utdt = U[i]

        # iteracion tipo Newton-Rapshon para equilibrio en t+Dt
        while (j < maxiter) && (Ferr > reltol)
            # calculate effective loads
            mᵢ = M * (a₀ * U[i] + a₂ * U′[i] + a₃ * U′′[i])
            cᵢ = C * (a₁ * U[i] + a₄ * U′[i] + a₅ * U′′[i])
            Feff = Fext(tvec[i] + Δt) + mᵢ + cᵢ - Fint(uktdt) - (a₀*M + a₁*C) * utdt
            Keff = KT(utdt) + a₀ * M + a₁ * C

            # solve for displacements
            du = Keff \ Feff
            utdt += du

            # TODO change stopping criterion
            Ferr = norm(Feff) / norm(Fext(t[k]+Δt))

            j += 1
        end

        if j == maxiter
            @warn("maximum number of iterations reached")
            break
        end

        U[i+1]   = utdt
        U′′[i+1] = a₀*(U[i+1] - U[i]) - a₂* U′[i] - a₃ * U′′[i]
        U′[i+1]  = U′[i] + a₆*U′′[i] + a₇*U′′[i+1]

        tvec[i+1] = tvec[i] + Δt
        i += 1
    end

    return _build_solution(alg, U, U′, U′′, NSTEPS)
end
