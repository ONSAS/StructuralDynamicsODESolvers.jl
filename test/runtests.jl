using Test, StructuralDynamicsODESolvers

# ---------------------------------------------------------
# The examples in this file can be found in
# [Chapter 9, in Finite Element Procedures, K-J Bathe].
# ---------------------------------------------------------

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
example_9_1_Bathe = SecondOrderAffineContinuousSystem(M, C, K, R)

@testset "Central difference method" begin
    # Ref. Example 9.1 pp. 773-774
    U₀ = zeros(2)
    U₀′ = zeros(2)
    alg = CentralDifference(Δt=0.28)
    prob = InitialValueProblem(example_9_1_Bathe, (U₀, U₀′))

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.02) < 5e-3
    @test abs(sol[13][2] - 2.60) < 5e-3
end

@testset "Houbolt method" begin
    # Ref. Example 9.2 pp. 776-777
    U₀ = zeros(2)
    U₀′ = zeros(2)
    alg = Houbolt(Δt=0.28)
    prob = InitialValueProblem(example_9_1_Bathe, (U₀, U₀′))

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.72) < 5e-3
    @test abs(sol[13][2] - 2.28) < 5e-3
end

@testset "Newmark method" begin
    # Ref. Example 9.3 pp. 778-779
    U₀ = zeros(2)
    U₀′ = zeros(2)
    alg = Newmark(Δt=0.28, α=0.25, δ=0.5)
    prob = InitialValueProblem(example_9_1_Bathe, (U₀, U₀′))

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.40) < 5e-3
    @test abs(sol[13][2] - 2.31) < 5e-3
end

@testset "Bathe method" begin
    # Ref. Example 9.4 pp. 781-782
    U₀ = zeros(2)
    U₀′ = zeros(2)
    alg = Bathe(Δt=0.28)
    prob = InitialValueProblem(example_9_1_Bathe, (U₀, U₀′))

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.28) < 5e-3
    @test abs(sol[13][2] - 2.40) < 5e-3
end

@testset "BackwardEuler method" begin

    C = [1. 0; 0 1.]
    K = [1. -1.; -1. 1.]
    M = zeros(2, 2)
    R = [0., 1.]
    heatTransferProblem = SecondOrderAffineContinuousSystem(M, C, K, R)
    alg = BackwardEuler(Δt=0.1)

    U₀ = zeros(2)
    prob = InitialValueProblem( heatTransferProblem, (U₀,U₀) )
    @test abs(0) < 1e-3
end
