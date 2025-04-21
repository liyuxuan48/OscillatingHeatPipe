@testset "Liquid arrays" begin

    L, Lmin, ϕ0 = 4.0, 0.001, 0.5
    closedornot = true
    X, dXdt, realratio = randomXp(L,Lmin,closedornot;chargeratio=ϕ0)

    l = mod.(map(u -> u[2],X) - map(u -> u[1],X),L)
    @test all(l .> 0)

    ϕ = sum(l)/L
    @test abs(ϕ - ϕ0) < 0.1
    @test ϕ ≈ realratio

end