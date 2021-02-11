example_9_1_Bathe = oscillator()
heatTransferProblem = heat_transfer()

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

    U₀ = zeros(2)
    for j=1:(nelem+1-2)
      U₀[j] =  sin( pi*lelem*j ) + 0.5 * sin( 3.0*pi*lelem*j )
    end
    testΔt = 0.0001
    # struct algorithm
    alg = BackwardEuler(Δt=testΔt)

    # struct ivp
    prob = InitialValueProblem( heatTransferProblem, (U₀,U₀) )

    α = kco / ( rho * csh )

    sol = solve(prob, alg, NSTEPS=12) |> displacements

    j=1
    analyticVal = exp(-pi^2 * α * 12*testΔt ) * sin( pi*lelem*j ) + 0.5 * exp(-(3*pi)^2 * α * 12*testΔt ) * sin( 3.0*pi*lelem*j )

    # relative error at node 1 verification
    @test ( abs( sol[13][1] - analyticVal ) / abs( analyticVal ) ) < 5e-3
end
