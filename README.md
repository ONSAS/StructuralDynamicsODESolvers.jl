# StructuralDynamicsODESolvers.jl

[![Build Status](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/workflows/CI/badge.svg)](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/actions?query=workflow%3ACI)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://onsas.github.io/StructuralDynamicsODESolvers.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/blob/master/LICENSE)
[![Join the chat at https://gitter.im/ONSAS_/community](https://badges.gitter.im/ONSAS_/community.svg)](https://gitter.im/ONSAS_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)


This package contains pure Julia implementations of ordinary differential equations (ODE) solvers for
structural dynamics problems.

## Features

The following solvers for linear dynamic equations are available:

- Bathe (equal size sub-steps) [BAT07]
- Central difference
- Houbolt [HOU50]
- Newmark [NEW509]
- Backward Euler (for first order systems)

## Example

The following example is explained in [this notebook](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/blob/master/examples/massDashpotSpringExample.ipynb).

For further examples see the [Example](https://nbviewer.jupyter.org/github/ONSAS/StructuralDynamicsODESolvers.jl/blob/gh-pages/dev/models/massDashpotSpring.ipynb) section of the documentation.

```julia
using StructuralDynamicsODESolvers

k  = 2 ; m  = .5 ;  c = .1 
u0 = 1 ; v0 = 0 

alg = Bathe(Δt = 0.1)

M = m*ones(1, 1)
C = c*ones(1, 1)
K = k*ones(1, 1)
R = zeros(1)

sys = SecondOrderAffineContinuousSystem(M, C, K, R)

U₀ = u0 * ones(1) ; V₀ = v0 * ones(1)

prob = InitialValueProblem(sys, (U₀, V₀))

sol = solve(prob, alg, NSTEPS=300);
```

```julia
using Plots

plot(sol, vars=(0, 1), xlab="time", ylab="displacement")
```

## Related libraries

This package has been created for research purposes. If you are new to numerically solving differential equations in Julia, we strongly suggest that you use the [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/) suite.

## References


- [BAT07] Bathe, Klaus-Jürgen. "[Conserving energy and momentum in nonlinear dynamics: a simple implicit time integration scheme.](https://www.sciencedirect.com/science/article/abs/pii/S0045794906003099)" Computers & structures 85.7-8 (2007): 437-445.
- [NEW59] Newmark, Nathan M. "[A method of computation for structural dynamics.](https://cedb.asce.org/CEDBsearch/record.jsp?dockey=0011858)" Journal of the engineering mechanics division 85.3 (1959): 67-94.
- [HOU50] Houbolt, John C. "[A recurrence matrix solution for the dynamic response of elastic aircraft.](https://arc.aiaa.org/doi/10.2514/8.1722)" Journal of the Aeronautical Sciences 17.9 (1950): 540-550.
