```@meta
EditURL = "<unknown>/examples/massDashpotSpring.jl"
```

Mass dashpot spring

[![](https://img.shields.io/badge/show-nbviewer-579ACA.svg)](<unknown>/models/massDashpotSpring.ipynb)

```@example massDashpotSpring
# Free oscillations

using StructuralDynamicsODESolvers, Plots

k  = 2 ; m  = .5 ;  c = .1 ;
u0 = 1 ; v0 = 0 ;

M = m * ones(1, 1)
C = c * ones(1, 1)
K = k * ones(1, 1)
R = zeros(1)

sys = SecondOrderAffineContinuousSystem(M, C, K, R)

U₀ = u0 * ones(1)
V₀ = v0 * ones(1)

ivp_free = InitialValueProblem(sys, (U₀, V₀))

alg = Bathe(Δt = 0.1)
sol = solve(ivp_free, alg, NSTEPS=100)
```

```@example massDashpotSpring
plot(times(sol), displacements(sol, 1))

# Forced oscillations

NSTEPS = 100
Δt = 0.1
ωf = k/(2m)
R = [[0.1 * sin(ωf * Δt * (i-1))] for i in 1:NSTEPS+1]

X = nothing # state constraints are ignored
B = ones(1, 1)
sys = SecondOrderConstrainedLinearControlContinuousSystem(M, C, K, B, X, R)

U₀ = u0 * ones(1)
V₀ = v0 * ones(1)

ivp_forced = InitialValueProblem(sys, (U₀, V₀))

alg = Bathe(Δt = 0.1)
sol = solve(ivp_forced, alg, NSTEPS=NSTEPS)
```

```@example massDashpotSpring
plot(times(sol), displacements(sol, 1))
```

