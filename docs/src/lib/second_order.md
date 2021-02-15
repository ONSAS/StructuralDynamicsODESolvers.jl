```@meta
DocTestSetup = :(using StructuralDynamicsODESolvers)
CurrentModule = StructuralDynamicsODESolvers
```

# Second-order problems

This section includes direct integration methods for linear dynamic equations
of the form

```math
    Mx''(t) + Cx'(t) + Kx(t) = F(t)
```
In the context of structural dynamics problems, $M$ is the mass matrix, $C$ is
the viscous damping matrix, $K$ is the stiffness matrix, $F$ is the vector of applied
forces, and $x(t)$, $x'(t)$ and $x''(t)$ are the displacement, velocity and
acceleration vectors, respectively.

The following algorithms are available:

- Central difference method.
- Houbolt's method.
- Newmark's method.
- Bathe integration method with equal-size substeps ($\gamma = 0.5$).

The theoretical description of such methods can be found in Chapter 9, [[BATHE]](@ref);
see also the references appearing in each docstring.

## Central difference

```@docs
CentralDifference
```

## Houbolt

```@docs
Houbolt
```

## Newmark

```@docs
Linear
Newmark
Trapezoidal
```

## Bathe

```@docs
Bathe
```
