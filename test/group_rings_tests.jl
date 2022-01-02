@testset "cyclic_group" begin
    n = 3
    Cₙ = cyclic_group(n)
    a, = Groups.gens(Cₙ)

    @test a^n == one(Cₙ)
    @test all(a^i * a^j == a^(i + j) for i in 2:n for j in 2:n)
end

@testset "group_ring" begin
    n = 15
    Cₙ = cyclic_group(n)
    a, = Groups.gens(Cₙ)

    RCₙ = LowCohomologySOS.group_ring(Cₙ, n)
    b = RCₙ.basis
    mstr = RCₙ.mstructure
    @test RCₙ.object == Cₙ
    @test Set(collect(b)) == Set([[a]; [a^i for i in 2:n]])
    @test mstr[getindex(b, a), getindex(b, a)] == getindex(b, a^2)
    @test all(
        mstr[getindex(b, a), getindex(b, a^j)] == getindex(b, a^(j + 1)) for
        j in 2:n
    )
    @test all(
        mstr[getindex(b, a^i), getindex(b, a^j)] == getindex(b, a^(i + j)) for
        i in 2:n for j in 2:n
    )

    RCₙ_star = LowCohomologySOS.group_ring(Cₙ, n, true)
    b_star = RCₙ_star.basis
    mstr_star = RCₙ_star.mstructure
    @test RCₙ_star.object == Cₙ
    @test Set(collect(b_star)) == Set([[a]; [a^i for i in 2:n]])
    @test mstr_star[getindex(b_star, a), getindex(b_star, a)] ==
          getindex(b_star, one(Cₙ))
    @test all(
        mstr_star[getindex(b_star, a), getindex(b_star, a^j)] ==
        getindex(b_star, a^(-1 + j)) for j in 2:n if -1 + j != 1
    )
    @test all(
        mstr_star[getindex(b_star, a^i), getindex(b_star, a^j)] ==
        getindex(b_star, a^(-i + j)) for i in 2:n for
        j in 2:n if -i + j != 1 && i - j != 1
    )

    A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    x, y, z = Groups.gens(F)
    half_basis = [x, x * y * z]

    RF = LowCohomologySOS.group_ring(F, half_basis)
    b_RF = RF.basis
    mstr_RF = RF.mstructure
    @test RF.object == F
    @test Set(collect(b_RF)) == Set([
        one(F),
        x,
        x^(-1),
        x * y * z,
        (x * y * z)^(-1),
        x^2,
        x^(-2),
        (x * y * z)^2,
        (x * y * z)^(-2),
        x^2 * y * z,
        x * (x * y * z)^(-1),
        y * z,
        x^(-1) * (x * y * z)^(-1),
        x * y * z * x,
        x * y * z * x^(-1),
        (x * y * z)^(-1) * x,
        (x * y * z)^(-1) * x^(-1),
    ])
    @test mstr_RF[getindex(b_RF, x), getindex(b_RF, x * y * z)] ==
          getindex(b_RF, x^2 * y * z)
    @test mstr_RF[getindex(b_RF, x * y * z), getindex(b_RF, x * y * z)] ==
          getindex(b_RF, (x * y * z)^2)
    @test mstr_RF[getindex(b_RF, (x * y * z)^(-1)), getindex(b_RF, x^(-1))] ==
          getindex(b_RF, (x * y * z)^(-1) * x^(-1))

    RF_star = LowCohomologySOS.group_ring(F, half_basis, true)
    b_RF_star = RF_star.basis
    mstr_RF_star = RF_star.mstructure
    @test RF_star.object == F
    @test Set(collect(b_RF_star)) == Set(collect(b_RF))
    @test mstr_RF_star[
        getindex(b_RF_star, x),
        getindex(b_RF_star, x * y * z),
    ] == getindex(b_RF_star, y * z)
    @test mstr_RF_star[
        getindex(b_RF_star, x * y * z),
        getindex(b_RF_star, x * y * z),
    ] == getindex(b_RF_star, one(F))
    @test mstr_RF_star[
        getindex(b_RF_star, (x * y * z)^(-1)),
        getindex(b_RF_star, x^(-1)),
    ] == getindex(b_RF_star, x * y * z * x^(-1))
end
