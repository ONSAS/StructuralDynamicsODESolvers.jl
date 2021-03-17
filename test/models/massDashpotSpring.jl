using StructuralDynamicsODESolvers, Plots, LinearAlgebra

k  = 2 ; m  = .5 ;  c = 0 ;
u0 = 1 ; v0 = 0 ;

M = m * ones(1, 1)
C = c * ones(1, 1)
K = k * ones(1, 1)
R = zeros(1)

sys = SecondOrderAffineContinuousSystem(M, C, K, R)

U₀ = u0 * ones(1); V₀ = v0 * ones(1);

ivp_free = InitialValueProblem(sys, (U₀, V₀))

NSTEPS = 1000 ;
Δt = 0.01 ;

alg = Bathe(Δt = Δt )
sol = solve(ivp_free, alg, NSTEPS=NSTEPS);

plot(sol, vars=(0, 1))

ωN = k/m
ωf = ωN * 2
Af = 10.0
R  = [ [ Af * sin(ωf * Δt * (i-1) ) ] for i in 1:NSTEPS+1];

X   = nothing # state constraints are ignored
B   = ones(1, 1)
sys = SecondOrderConstrainedLinearControlContinuousSystem(M, C, K, B, X, R)

ivp_forced_secOrder = InitialValueProblem(sys, (U₀, V₀))

alg = Bathe(Δt = Δt )
sol_secOrder = solve(ivp_forced_secOrder, alg, NSTEPS=NSTEPS);

#The new vector of variables is

K = [     0 1     0 0 ;
      -ωN^2 0     1 0 ;
          0 0     0 1 ;
          0 0 -ωf^2 0 ] ;

C = -Diagonal(ones(4))
M = zeros(4,4)
R = zeros(4)

sys = SecondOrderAffineContinuousSystem(M, C, K, R)

U₀ = [u0; v0; 0; ωf*Af ] ;

ivp_forced_firOrder = InitialValueProblem(sys, (U₀, U₀) )

alg = BackwardEuler(Δt = Δt )
sol_firOrder = solve(ivp_forced_firOrder, alg, NSTEPS=NSTEPS);

plot(sol_secOrder, vars=(0, 1), xlab="time" )
plot!(sol_firOrder, vars=(0, 1), xlab="time" )

