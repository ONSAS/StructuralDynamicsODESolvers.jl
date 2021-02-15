```@meta
DocTestSetup = :(using StructuralDynamicsODESolvers)
CurrentModule = StructuralDynamicsODESolvers
```

# API Reference

## Solve interface

```@docs
Solution
step_size
solve
displacements
velocities
accelerations
```

## Type hierarchy & dispatch

Each integration algorithm, called *solver*, is implemented by a Julia struct
which holds the algorithm's options (such as the step-size), and such struct is
a subtype of `AbstractSolver`. New solvers should implement a method with signature

```julia
_solve(::SolverType, ::InitialValueProblem{...}, args...; kwargs...)`
```
The return type of a solver must be any concrete subtype that implements the
`AbstractSolution` interface, e.g. the `Solution` type.

```@docs
StructuralDynamicsODESolvers.AbstractSolution
StructuralDynamicsODESolvers.AbstractSolver
```
