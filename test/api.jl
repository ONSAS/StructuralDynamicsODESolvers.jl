@testset "Solve getter functions" begin
    ivp = harmonic_oscillator_free()

    alg = Bathe(; Δt=0.1)
    sol = solve(ivp, alg; NSTEPS=100)

    @test dim(sol) == 1

    U = [ui[1] for ui in sol.U]
    @test U == displacements(sol, 1)
    @test_throws ArgumentError displacements(sol, 2)

    V = [vi[1] for vi in sol.U′]
    @test V == velocities(sol, 1)
    @test_throws ArgumentError velocities(sol, 2)

    A = [ai[1] for ai in sol.U′′]
    @test A == accelerations(sol, 1)
    @test_throws ArgumentError accelerations(sol, 2)
end

@testset "Solver options" begin
    ivp = harmonic_oscillator_free()

    alg = Bathe(; Δt=0.1)

    sol = solve(ivp, alg; NSTEPS=100)
    @test times(sol)[end] == 10.0

    sol = solve(ivp, alg; T=10.0)
    @test times(sol)[end] == 10.0

    sol = solve(ivp, alg; tspan=(0.0, 10.0))
    @test times(sol)[end] == 10.0

    sol = solve(ivp, alg; finalTime=10.0)
    @test times(sol)[end] == 10.0

    # test unknown argument name
    @test_throws ArgumentError solve(ivp, alg, final_time=10.0)

    # starting with shifted interval is not implemented
    @test_throws AssertionError solve(ivp, alg, tspan=(1.0, 2.0))
end
