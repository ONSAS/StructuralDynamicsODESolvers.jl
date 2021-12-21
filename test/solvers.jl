@testset "Central difference method" begin
    # Ref. Example 9.1 pp. 773-774
    prob, U = example_9_1_Bathe()
    alg = CentralDifference(Δt=0.28)

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.02) < 5e-3
    @test abs(sol[13][2] - 2.60) < 5e-3
end

@testset "Houbolt method" begin
    # Ref. Example 9.2 pp. 776-777
    prob, U = example_9_1_Bathe()
    alg = Houbolt(Δt=0.28)

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.72) < 5e-3
    @test abs(sol[13][2] - 2.28) < 5e-3
end

@testset "Newmark method" begin
    # Ref. Example 9.3 pp. 778-779
    prob, U = example_9_1_Bathe()
    alg = Newmark(Δt=0.28, α=0.25, δ=0.5)

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.40) < 5e-3
    @test abs(sol[13][2] - 2.31) < 5e-3
end

@testset "Bathe method" begin
    # Ref. Example 9.4 pp. 781-782
    prob, U = example_9_1_Bathe()
    alg = Bathe(Δt=0.28)

    # solution
    sol = solve(prob, alg, NSTEPS=12) |> displacements

    @test abs(sol[13][1] - 1.28) < 5e-3
    @test abs(sol[13][2] - 2.40) < 5e-3
end

@testset "Bathe method with forcing term" begin
    ivp, NSTEPS, Δt = harmonic_oscillator_forced()

    alg = Bathe(Δt = Δt)
    sol = solve(ivp, alg, NSTEPS=NSTEPS)

    # TODO add analytic solution for test
end

@testset "BackwardEuler method" begin
    prob, U = heat_transfer()

    testΔt = 0.0001
    alg = BackwardEuler(Δt=testΔt)

    sol = solve(prob, alg, NSTEPS=12) |> displacements

    # relative error at node 1 verification
    analyticVal = U(testΔt)
    @test (abs(sol[13][1] - analyticVal) / abs(analyticVal)) < 5e-3
end

@testset "Nonlinear Newmark method: Pendulum" begin

    # load model
    prob = pendulum_nonlinear()

    # solve using Trapezoidal scheme
    α = 1/4
    δ = 1/2
    Δt = 0.01
    alg = Newmark(Δt=Δt, α=α, δ=δ)

    sol = solve(prob, alg, T=3.0)

    # FIXME review tests
    @test sol.U[end] ≈ [1.9113595927872666, 1.410419996390988]
    @test sol.U′[end] ≈ [1.0007721533056833, 3.2500475651804632]
    @test sol.U′′[end] ≈ [-8.286971439792033, -7.2537920159991955]
end

@testset "Nonlinear Newmark method: Von Mises Truss" begin

    # load model
    prob = von_mises_truss()

    # algorithm parameters
    α = 1/4
    δ = 1/2
    Δt = 0.0025
    alg = Newmark(Δt=Δt, α=α, δ=δ)

    # solve using Trapezoidal scheme
    sol = solve(prob, alg, T=2.0)

    # FIXME update tests
    sol.U[end] ≈ [0.001888137400895257, -0.17208311613153432]
    sol.U′[end] ≈ [-0.07080458802627602, 0.14912664021827504]
    sol.U′′[end] ≈ [6.0851727632872485, -6.18000018114007]
end

@testset "Nonlinear Central difference method: Von Mises Truss" begin

    # load model
    prob = von_mises_truss()

    # algorithm parameters
    Δt = 0.000025
    alg = CentralDifference(Δt=Δt)

    # solve
    sol = solve(prob, alg, T=2.0)

    # FIXME update tests
end
