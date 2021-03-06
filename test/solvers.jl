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
    @test ( abs( sol[13][1] - analyticVal ) / abs( analyticVal ) ) < 5e-3
end
