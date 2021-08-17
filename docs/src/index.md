# StructuralDynamicsODESolvers.jl

[![Build Status](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/workflows/CI/badge.svg)](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/actions?query=workflow%3ACI)
[![Documentation](https://img.shields.io/badge/docs-latest-blue.svg)](https://onsas.github.io/StructuralDynamicsODESolvers.jl/dev/)
[![license](https://img.shields.io/github/license/mashape/apistatus.svg?maxAge=2592000)](https://github.com/ONSAS/StructuralDynamicsODESolvers.jl/blob/master/LICENSE)
[![Join the chat at https://gitter.im/ONSAS_/community](https://badges.gitter.im/ONSAS_/community.svg)](https://gitter.im/ONSAS_/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

This package contains pure Julia implementations of ordinary differential equations (ODE)
solvers for structural dynamics problems.

## Features

The following algorithms for linear dynamic equations are implemented:

For **second order**  problems:

- Bathe (equal size sub-steps) [[BAT07]](@ref)
- Central difference
- Houbolt [[HOU50]](@ref)
- Newmark [[NEW59]](@ref)

For **first order** problems:

- Backward Euler

## Related libraries

This package has been created for research purposes. If you are new to numerically solving differential equations in Julia, we strongly suggest that you use the [DifferentialEquations.jl](https://diffeq.sciml.ai/dev/) suite.

## References

See the [References](@ref refs_page) section.

## Contents

```@contents
Pages = [
    "lib/first_order.md",
    "lib/second_order.md",
    "lib/example.md",
    "lib/api.md"
]
Depth = 2
```
