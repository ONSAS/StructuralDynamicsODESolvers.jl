# ---------------------------------------------------------
# Simple harmonic oscillator without forcing term
# ---------------------------------------------------------
function harmonic_oscillator_free()
    k  = 2 ; m  = .5 ;  c = .1 ;
    u0 = 1 ; v0 = 0 ;

    M = m * ones(1, 1)
    C = c * ones(1, 1)
    K = k * ones(1, 1)
    R = zeros(1)

    sys = SecondOrderAffineContinuousSystem(M, C, K, R)

    U₀ = u0 * ones(1)
    V₀ = v0 * ones(1)

    return InitialValueProblem(sys, (U₀, V₀))
end

# ---------------------------------------------------------
# Simple harmonic oscillator with forcing term
# ---------------------------------------------------------
function harmonic_oscillator_forced()
    k  = 2 ; m  = .5 ;  c = .1 ;
    u0 = 1 ; v0 = 0 ;

    M = m * ones(1, 1)
    C = c * ones(1, 1)
    K = k * ones(1, 1)

    NSTEPS = 100
    Δt = 0.1
    ωf = k/(2m)
    R = [[0.1 * sin(ωf * Δt * (i-1))] for i in 1:NSTEPS+1]

    X = nothing # Universe(1) ignores state constraints
    B = ones(1, 1)
    sys = SecondOrderConstrainedLinearControlContinuousSystem(M, C, K, B, X, R)

    U₀ = u0 * ones(1)
    V₀ = v0 * ones(1)

    return InitialValueProblem(sys, (U₀, V₀)), NSTEPS, Δt
end

# ---------------------------------------------------------
# This example can be found in
# [Chapter 9, in Finite Element Procedures, K-J Bathe].
# ---------------------------------------------------------
function example_9_1_Bathe(; U₀ = zeros(2),   # initial displacements
                             U₀′ = zeros(2)) # initial velocities

    # analytic solution (cf. Example 9.7)
    A = [1/√3  (1/2)*√(2/3);
         1/√3      -√(2/3)]
    x₁(t) = (5 / √3) * (1 - cos(t*√2))
    x₂(t) = (2 * √(2/3)) * (-1 + cos(t*√5))
    U(t) = A * [x₁(t), x₂(t)]

    U12 = U(12*0.28) # ≈ [1.157, 2.489]

    # problem formulation
    M = [2 0; 0 1.]
    K = [6 -2; -2 4.]
    C = zeros(2, 2)
    R = [0, 10.]

    sys = SecondOrderAffineContinuousSystem(M, C, K, R)
    return InitialValueProblem(sys, (U₀, U₀′)), U
end

# ---------------------------------------------------------
# Heat transfer problem
# ---------------------------------------------------------
function heat_transfer(; nelem=3)
    rho = 1.
    csh = 1.
    kco = 4.
    lelem = 1.0 / nelem
    Area  = .25
    C = rho * csh * lelem * Area * ( [0.5+0.5+0.0 0; 0  0.0+0.5+0.5] )
    K = kco *       lelem * Area * ( [1.0+1.0 -1.; -1. 1.0+1.0] * 1.0 / lelem^2 )
    M = zeros(2, 2)
    R = zeros(2)
    sys = SecondOrderAffineContinuousSystem(M, C, K, R)

    # build initial state
    U₀ = zeros(2)
    for j=1:(nelem+1-2)
      U₀[j] =  sin( pi*lelem*j ) + 0.5 * sin( 3.0*pi*lelem*j )
    end

    # analytic solution (node 1)
    α = kco / ( rho * csh )
    j = 1
    U(t) = exp(-pi^2 * α * 12*t ) * sin( pi*lelem*j ) + 0.5 * exp(-(3*pi)^2 * α * 12*t ) * sin( 3.0*pi*lelem*j )

    # here the velocity condition U₀′ is ignored
    return InitialValueProblem(sys, (U₀, U₀)), U
end
