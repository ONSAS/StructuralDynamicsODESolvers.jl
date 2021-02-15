```@meta
DocTestSetup = :(using StructuralDynamicsODESolvers)
CurrentModule = StructuralDynamicsODESolvers
```

# First-order problems

This section includes direct integration methods for linear dynamic equations
of the form:

```math
    Mx'(t) + Kx(t) = F(t)
```
In the context of heat transfer problems,  $M$ is the capacity matrix, $K$ is
the conductivity matrix, $F(t)$ is the heat supply vector, $x(t)$ is the temperature
vector, and $x'(t)$ is the time derivative of $x(t)$.

The following algorithms are available:

- Backward (implicit) Euler

## Backward Euler

```@docs
BackwardEuler
```
