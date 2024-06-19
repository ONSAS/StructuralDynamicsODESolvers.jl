# ---------------------------------------------------------
# Simple harmonic oscillator without forcing term
# ---------------------------------------------------------
function harmonic_oscillator_free()
    k = 2
    m = 0.5
    c = 0.1
    u0 = 1
    v0 = 0

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
    k = 2
    m = 0.5
    c = 0.1
    u0 = 1
    v0 = 0

    M = m * ones(1, 1)
    C = c * ones(1, 1)
    K = k * ones(1, 1)

    NSTEPS = 100
    Δt = 0.1
    ωf = k / (2m)
    R = [[0.1 * sin(ωf * Δt * (i - 1))] for i in 1:(NSTEPS + 1)]

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
function example_9_1_Bathe(;
                           U₀=zeros(2),   # initial displacements
                           U₀′=zeros(2),) # initial velocities

    # analytic solution (cf. Example 9.7)
    A = [1/√3 (1 / 2)*√(2 / 3)
         1/√3 -√(2 / 3)]
    x₁(t) = (5 / √3) * (1 - cos(t * √2))
    x₂(t) = (2 * √(2 / 3)) * (-1 + cos(t * √5))
    U(t) = A * [x₁(t), x₂(t)]

    U12 = U(12 * 0.28) # ≈ [1.157, 2.489]

    # problem formulation
    M = [2 0; 0 1.0]
    K = [6 -2; -2 4.0]
    C = zeros(2, 2)
    R = [0, 10.0]

    sys = SecondOrderAffineContinuousSystem(M, C, K, R)
    return InitialValueProblem(sys, (U₀, U₀′)), U
end

# ---------------------------------------------------------
# Heat transfer problem
# ---------------------------------------------------------
function heat_transfer(; nelem=3)
    rho = 1.0
    csh = 1.0
    kco = 4.0
    lelem = 1.0 / nelem
    Area = 0.25
    C = rho * csh * lelem * Area * ([0.5+0.5+0.0 0; 0 0.0+0.5+0.5])
    K = kco * lelem * Area * ([1.0+1.0 -1.0; -1.0 1.0+1.0] * 1.0 / lelem^2)
    M = zeros(2, 2)
    R = zeros(2)
    sys = SecondOrderAffineContinuousSystem(M, C, K, R)

    # build initial state
    U₀ = zeros(2)
    for j in 1:(nelem + 1 - 2)
        U₀[j] = sin(pi * lelem * j) + 0.5 * sin(3.0 * pi * lelem * j)
    end

    # analytic solution (node 1)
    α = kco / (rho * csh)
    j = 1
    function U(t)
        exp(-pi^2 * α * 12 * t) * sin(pi * lelem * j) +
        0.5 * exp(-(3 * pi)^2 * α * 12 * t) * sin(3.0 * pi * lelem * j)
    end

    # here the velocity condition U₀′ is ignored
    return InitialValueProblem(sys, (U₀, U₀)), U
end

# ---------------------------------------------------------
# Non-linear pendulum model
# ---------------------------------------------------------
function pendulum_nonlinear()
    # ========================
    # Model definition
    # ========================

    # bar length
    L = 2.0 # [m]

    # cross section area
    ϕ = 0.01 # [m]
    A = (π * ϕ^2 / 4) # [m^2]

    # Young's modulus (Steel)
    E = 210e9 # [Pa]

    # mass (lumped)
    m = 214 # [kg]

    # damping
    c = 0 # [kg/s]

    # gravity
    g = 9.81 # [m/s^2]

    # ============================
    # Finite element model
    # ============================

    # u = [u_1, u_2]^T is the position of the mass, with the origin as indicated:
    #                     _______
    #                        |
    #                        |
    #                        |
    #                        o    << origin

    # internal forces field
    daux(u) = u[1]^2 - 2 * L * u[2] + u[2]^2
    Fint(u) = E * A * L * daux(u) / 2 / L^4 * [u[1], -(L - u[2])]

    # external forces field
    Fext(t) = [0, -m * g] # [N]

    # lumped mass matrix
    M = [m 0; 0.0 m]

    # damping matrix
    C = [c 0; 0.0 c]

    # ===================================
    # Initial conditions for simulation
    # ===================================

    u0 = [L, L]
    v0 = [0, 0.0]

    sys = SecondOrderContinuousSystem(M, C, Fint, Fext)
    prob = InitialValueProblem(sys, (u0, v0))

    return prob
end

# ---------------------------------------------------------
# Von Mises truss (cercha de Von Mises)
# ---------------------------------------------------------
function von_mises_truss()
    # structural parameters
    ρ = 7850               # [kg/m^3] (steel)
    Lx = 0.374 / 2
    Lz = sqrt(0.205^2 - Lx^2) # [m]
    l0 = sqrt(Lx^2 + Lz^2) # [m]
    Lc = 0.240              # [m]
    Ic = 0.0254 * 0.0032^3 / 12  # [m^4]
    Ac = 0.0254 * 0.0032       # [m^2]
    E = 200e9             # [Pa] (steel)
    kc = 3 * E * Ic / Lc^3       # [N/m]
    mb = l0 * Ac * ρ           # [kg]
    m = 1.4               # kg (pandeo incipiente en 1.4)
    c = 2                 # [kg/s] (amortiguamiento por friccion juntas y arrastre pesa)
    g = 9.81              # [m/s^2]

    # internal forces vector
    function Fint(u)
        E * Ac * l0 * (u[1]^2 + u[2]^2 - 2 * Lx * u[1] + 2 * Lz * u[2]) / 2 / l0^4 *
        [-Lx + u[1], Lz + u[2]] + [kc * u[1], 0]
    end

    # external forces vector
    Fext(t) = [0, -(m + mb) / 2 * g] # [N]

    # mass matrix (lumped)
    M = [mb 0; 0 (mb + m)/2]

    # damping matrix
    C = [c/10 0; 0 c]

    # initial conditions for simulation
    u0 = [0, 0.0]
    v0 = [0, 0.0]

    sys = SecondOrderContinuousSystem(M, C, Fint, Fext)
    prob = InitialValueProblem(sys, (u0, v0))
end
