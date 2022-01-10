@testset "d₀" begin
    A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    RF = LowCohomologySOS.group_ring(F, 1)
    x, y, z = Groups.gens(F)
    @test LowCohomologySOS.d₀(RF, [x, y, z]) ==
          reshape([RF(x) - one(RF), RF(y) - one(RF), RF(z) - one(RF)], 3, 1)

    G = FPGroup(F, [x * y => y * x, y * z => z * y, z * x => x * z])
    RG = LowCohomologySOS.group_ring(G, 1)
    xx, yy, zz = Groups.gens(G)
    @test LowCohomologySOS.d₀(RG, [xx, yy, zz]) ==
          reshape([RG(xx) - one(RG), RG(yy) - one(RG), RG(zz) - one(RG)], 3, 1)
end

@testset "embed" begin
    A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F = FreeGroup(A)
    RF = LowCohomologySOS.group_ring(F, 1)
    x, y = Groups.gens(F)
    G = FPGroup(F, [x * y => y * x])
    RG = LowCohomologySOS.group_ring(G, 1)
    xx, yy = Groups.gens(G)

    function zero_hom(u::FPGroupElement)
        return one(G)
    end
    @test LowCohomologySOS.embed(zero_hom, zero(RF), RG) == zero(RG)
    @test LowCohomologySOS.embed(zero_hom, RF(x), RG) == one(RG)
    @test LowCohomologySOS.embed(zero_hom, RF(y), RG) == one(RG)
    @test LowCohomologySOS.embed(zero_hom, RF(x * y) + RF(y^2), RG) ==
          2 * one(RG)

    @test LowCohomologySOS.embed(x -> one(G), zero(RF), RG) == zero(RG)
    @test LowCohomologySOS.embed(x -> one(G), RF(x), RG) == one(RG)
    @test LowCohomologySOS.embed(x -> one(G), RF(y), RG) == one(RG)
    @test LowCohomologySOS.embed(x -> one(G), RF(x * y) + RF(y^2), RG) ==
          2 * one(RG)

    function quotient_hom(u::FPGroupElement)
        result = one(G)
        for i in 1:length(word(u))
            if F(word(u)[i:i]) == x
                result *= xx
            elseif F(word(u)[i:i]) == inv(x)
                result *= inv(xx)
            elseif F(word(u)[i:i]) == y
                result *= yy
            elseif F(word(u)[i:i]) == inv(y)
                result *= inv(yy)
            end
        end
        return result
    end
    @test LowCohomologySOS.embed(quotient_hom, zero(RF), RG) == zero(RG)
    @test LowCohomologySOS.embed(quotient_hom, RF(x), RG) == RG(xx)
    @test LowCohomologySOS.embed(quotient_hom, RF(y), RG) == RG(yy)
    @test LowCohomologySOS.embed(quotient_hom, RF(x * y) + RF(y^2), RG) ==
          RG(xx * yy) + RG(yy^2)

    A2 = Alphabet(
        [
            :e12,
            :E12,
            :e21,
            :E21,
            :e13,
            :E13,
            :e31,
            :E31,
            :e23,
            :E23,
            :e32,
            :E32,
        ],
        [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11],
    )
    F_sl_3_z = FreeGroup(A2)
    RF_sl_3_z = LowCohomologySOS.group_ring(F_sl_3_z, 1)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F_sl_3_z)
    N = 3

    SL₃Ƶ = AbstractAlgebra.MatrixAlgebra(AbstractAlgebra.zz, N)
    StarAlgebras.star(A::AbstractAlgebra.Generic.MatAlgElem) = inv(A)

    E(SL₃Ƶ, i, j) = (e_ij = one(SL₃Ƶ);
    e_ij[i, j] = 1;
    e_ij)
    S = [E(SL₃Ƶ, i, j) for i in 1:N for j in 1:N if i ≠ j]
    half_basis, sizes = Groups.wlmetric_ball(unique([S; inv.(S)]), radius = 1)
    RSL₃Ƶ = LowCohomologySOS.group_ring(SL₃Ƶ, half_basis, true)
    function h(u::FPGroupElement)
        result = one(SL₃Ƶ)
        for i in 1:length(word(u))
            if F_sl_3_z(word(u)[i:i]) == e12
                result *= S[1]
            elseif F_sl_3_z(word(u)[i:i]) == inv(e12)
                result *= inv(S[1])
            elseif F_sl_3_z(word(u)[i:i]) == e21
                result *= S[3]
            elseif F_sl_3_z(word(u)[i:i]) == inv(e21)
                result *= inv(S[3])
            elseif F_sl_3_z(word(u)[i:i]) == e13
                result *= S[2]
            elseif F_sl_3_z(word(u)[i:i]) == inv(e13)
                result *= inv(S[2])
            elseif F_sl_3_z(word(u)[i:i]) == e31
                result *= S[5]
            elseif F_sl_3_z(word(u)[i:i]) == inv(e31)
                result *= inv(S[5])
            elseif F_sl_3_z(word(u)[i:i]) == e23
                result *= S[4]
            elseif F_sl_3_z(word(u)[i:i]) == inv(e23)
                result *= inv(S[4])
            elseif F_sl_3_z(word(u)[i:i]) == e32
                result *= S[6]
            elseif F_sl_3_z(word(u)[i:i]) == inv(e32)
                result *= inv(S[6])
            end
        end
        return result
    end
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e12), RSL₃Ƶ) == RSL₃Ƶ(S[1])
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e21), RSL₃Ƶ) == RSL₃Ƶ(S[3])
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e13), RSL₃Ƶ) == RSL₃Ƶ(S[2])
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e31), RSL₃Ƶ) == RSL₃Ƶ(S[5])
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e23), RSL₃Ƶ) == RSL₃Ƶ(S[4])
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e32), RSL₃Ƶ) == RSL₃Ƶ(S[6])
    @test LowCohomologySOS.embed(h, RF_sl_3_z(e12 * e13^(-1)), RSL₃Ƶ) ==
          RSL₃Ƶ(S[1] * S[2]^(-1))
end

@testset "fox_derivative" begin
    A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F = FreeGroup(A)
    RF = LowCohomologySOS.group_ring(F, 2)

    @test LowCohomologySOS.fox_derivative(RF, one(F), 1) == zero(RF)
    @test LowCohomologySOS.fox_derivative(RF, one(F), 2) == zero(RF)

    x, y = Groups.gens(F)

    @test LowCohomologySOS.fox_derivative(RF, x, 1) == one(RF)
    @test LowCohomologySOS.fox_derivative(RF, x, 2) == zero(RF)
    @test LowCohomologySOS.fox_derivative(RF, y, 1) == zero(RF)
    @test LowCohomologySOS.fox_derivative(RF, y, 2) == one(RF)

    @test LowCohomologySOS.fox_derivative(RF, x^(-1), 1) == -RF(x^(-1))
    @test LowCohomologySOS.fox_derivative(RF, y^(-1), 2) == -RF(y^(-1))

    @test LowCohomologySOS.fox_derivative(RF, x * y, 1) == one(RF)
    @test LowCohomologySOS.fox_derivative(RF, x * y, 2) == RF(x)

    @test LowCohomologySOS.fox_derivative(RF, x * y * x * y, 1) ==
          one(RF) + RF(x * y)
    @test LowCohomologySOS.fox_derivative(RF, x * y * x * y, 2) ==
          RF(x) + RF(x * y * x)
end

@testset "jacobian_matrix" begin
    A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F = FreeGroup(A)
    x, y = Groups.gens(F)
    relations = [x, x * y, x * y * x^(-1)]
    J = LowCohomologySOS.jacobian_matrix(relations)
    RF_J = parent(rand(J))
    RF_J_proper = LowCohomologySOS.suitable_group_ring(relations)
    @test RF_J.object == RF_J_proper.object
    @test RF_J.basis == RF_J_proper.basis
    @test RF_J.mstructure == RF_J_proper.mstructure
    @test typeof(RF_J) == typeof(RF_J_proper)
    J_proper = reshape([zero(RF_J) for i in 1:6], 3, 2)
    J_proper[1, 1] = one(RF_J)
    J_proper[1, 2] = zero(RF_J)
    J_proper[2, 1] = one(RF_J)
    J_proper[2, 2] = RF_J(x)
    J_proper[3, 1] = one(RF_J) - RF_J(x * y * x^(-1))
    J_proper[3, 2] = RF_J(x)
    @test J == J_proper

    A2 = Alphabet(
        [
            :e12,
            :E12,
            :e21,
            :E21,
            :e13,
            :E13,
            :e31,
            :E31,
            :e23,
            :E23,
            :e32,
            :E32,
        ],
        [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11],
    )
    F2 = FreeGroup(A2)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F2)
    relations2 = [
        e12 * e13 * e12^(-1) * e13^(-1),
        e12 * e32 * e12^(-1) * e32^(-1),
        e13 * e23 * e13^(-1) * e23^(-1),
        e23 * e21 * e23^(-1) * e21^(-1),
        e21 * e31 * e21^(-1) * e31^(-1),
        e31 * e32 * e31^(-1) * e32^(-1),
        e12 * e23 * e12^(-1) * e23^(-1) * e13^(-1),
        e13 * e32 * e13^(-1) * e32^(-1) * e12^(-1),
        e21 * e13 * e21^(-1) * e13^(-1) * e23^(-1),
        e23 * e31 * e23^(-1) * e31^(-1) * e21^(-1),
        e31 * e12 * e31^(-1) * e12^(-1) * e32^(-1),
        e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1),
    ]
    J2 = LowCohomologySOS.jacobian_matrix(relations2)
    RF2_J = parent(rand(J2))
    J2_proper = reshape([zero(RF2_J) for i in 1:72], 12, 6)
    J2_proper[1, 1] = one(RF2_J) - RF2_J(e12 * e13 * e12^(-1))
    J2_proper[1, 2] = zero(RF2_J)
    J2_proper[1, 3] = RF2_J(e12) - RF2_J(e12 * e13 * e12^(-1) * e13^(-1))
    J2_proper[1, 4] = zero(RF2_J)
    J2_proper[1, 5] = zero(RF2_J)
    J2_proper[1, 6] = zero(RF2_J)
    J2_proper[2, 1] = one(RF2_J) - RF2_J(e12 * e32 * e12^(-1))
    J2_proper[2, 2] = zero(RF2_J)
    J2_proper[2, 3] = zero(RF2_J)
    J2_proper[2, 4] = zero(RF2_J)
    J2_proper[2, 5] = zero(RF2_J)
    J2_proper[2, 6] = RF2_J(e12) - RF2_J(e12 * e32 * e12^(-1) * e32^(-1))
    J2_proper[3, 1] = zero(RF2_J)
    J2_proper[3, 2] = zero(RF2_J)
    J2_proper[3, 3] = one(RF2_J) - RF2_J(e13 * e23 * e13^(-1))
    J2_proper[3, 4] = zero(RF2_J)
    J2_proper[3, 5] = RF2_J(e13) - RF2_J(e13 * e23 * e13^(-1) * e23^(-1))
    J2_proper[3, 6] = zero(RF2_J)
    J2_proper[4, 1] = zero(RF2_J)
    J2_proper[4, 2] = RF2_J(e23) - RF2_J(e23 * e21 * e23^(-1) * e21^(-1))
    J2_proper[4, 3] = zero(RF2_J)
    J2_proper[4, 4] = zero(RF2_J)
    J2_proper[4, 5] = one(RF2_J) - RF2_J(e23 * e21 * e23^(-1))
    J2_proper[4, 6] = zero(RF2_J)
    J2_proper[5, 1] = zero(RF2_J)
    J2_proper[5, 2] = one(RF2_J) - RF2_J(e21 * e31 * e21^(-1))
    J2_proper[5, 3] = zero(RF2_J)
    J2_proper[5, 4] = RF2_J(e21) - RF2_J(e21 * e31 * e21^(-1) * e31^(-1))
    J2_proper[5, 5] = zero(RF2_J)
    J2_proper[5, 6] = zero(RF2_J)
    J2_proper[6, 1] = zero(RF2_J)
    J2_proper[6, 2] = zero(RF2_J)
    J2_proper[6, 3] = zero(RF2_J)
    J2_proper[6, 4] = one(RF2_J) - RF2_J(e31 * e32 * e31^(-1))
    J2_proper[6, 5] = zero(RF2_J)
    J2_proper[6, 6] = RF2_J(e31) - RF2_J(e31 * e32 * e31^(-1) * e32^(-1))
    J2_proper[7, 1] = one(RF2_J) - RF2_J(e12 * e23 * e12^(-1))
    J2_proper[7, 2] = zero(RF2_J)
    J2_proper[7, 3] = -RF2_J(e12 * e23 * e12^(-1) * e23^(-1) * e13^(-1))
    J2_proper[7, 4] = zero(RF2_J)
    J2_proper[7, 5] = RF2_J(e12) - RF2_J(e12 * e23 * e12^(-1) * e23^(-1))
    J2_proper[7, 6] = zero(RF2_J)
    J2_proper[8, 1] = -RF2_J(e13 * e32 * e13^(-1) * e32^(-1) * e12^(-1))
    J2_proper[8, 2] = zero(RF2_J)
    J2_proper[8, 3] = one(RF2_J) - RF2_J(e13 * e32 * e13^(-1))
    J2_proper[8, 4] = zero(RF2_J)
    J2_proper[8, 5] = zero(RF2_J)
    J2_proper[8, 6] = RF2_J(e13) - RF2_J(e13 * e32 * e13^(-1) * e32^(-1))
    J2_proper[9, 1] = zero(RF2_J)
    J2_proper[9, 2] = one(RF2_J) - RF2_J(e21 * e13 * e21^(-1))
    J2_proper[9, 3] = RF2_J(e21) - RF2_J(e21 * e13 * e21^(-1) * e13^(-1))
    J2_proper[9, 4] = zero(RF2_J)
    J2_proper[9, 5] = -RF2_J(e21 * e13 * e21^(-1) * e13^(-1) * e23^(-1))
    J2_proper[9, 6] = zero(RF2_J)
    J2_proper[10, 1] = zero(RF2_J)
    J2_proper[10, 2] = -RF2_J(e23 * e31 * e23^(-1) * e31^(-1) * e21^(-1))
    J2_proper[10, 3] = zero(RF2_J)
    J2_proper[10, 4] = RF2_J(e23) - RF2_J(e23 * e31 * e23^(-1) * e31^(-1))
    J2_proper[10, 5] = one(RF2_J) - RF2_J(e23 * e31 * e23^(-1))
    J2_proper[10, 6] = zero(RF2_J)
    J2_proper[11, 1] = RF2_J(e31) - RF2_J(e31 * e12 * e31^(-1) * e12^(-1))
    J2_proper[11, 2] = zero(RF2_J)
    J2_proper[11, 3] = zero(RF2_J)
    J2_proper[11, 4] = one(RF2_J) - RF2_J(e31 * e12 * e31^(-1))
    J2_proper[11, 5] = zero(RF2_J)
    J2_proper[11, 6] = -RF2_J(e31 * e12 * e31^(-1) * e12^(-1) * e32^(-1))
    J2_proper[12, 1] = zero(RF2_J)
    J2_proper[12, 2] = RF2_J(e32) - RF2_J(e32 * e21 * e32^(-1) * e21^(-1))
    J2_proper[12, 3] = zero(RF2_J)
    J2_proper[12, 4] = -RF2_J(e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1))
    J2_proper[12, 5] = zero(RF2_J)
    J2_proper[12, 6] = one(RF2_J) - RF2_J(e32 * e21 * e32^(-1))
    @test J2 == J2_proper
end

@testset "print_matrix" begin
    M = reshape([i for i in 1:9], 3, 3)
    f(io, m::AbstractMatrix) = Base.print_matrix(io, m)
    @test Base.sprint(f, M) == " 1  4  7\n 2  5  8\n 3  6  9"

    A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F = FreeGroup(A)
    x, y = Groups.gens(F)
    RF = LowCohomologySOS.group_ring(F, 1)

    Mx = reshape([zero(RF) for i in 1:4], 2, 2)
    Mx[1, 1] = zero(RF)
    Mx[1, 2] = one(RF)
    Mx[2, 1] = RF(x)
    Mx[2, 2] = RF(y)

    @test Base.sprint(f, Mx) == " 0·(id)  1·(id)\n 1·x     1·y"

    Mxx = copy(Mx)
    Mxx[1, 1] = one(RF) + RF(x * y^(-1))
    Mxx[1, 2] = 2 * RF(y * x)
    Mxx[2, 1] = RF(y^(-1))
    Mxx[2, 2] = RF(x^(-1))
    @test Base.sprint(f, Mxx) == " 1·(id) +1·x*Y  2·y*x\n 1·Y            1·X"
end

@testset "suitable_group_ring" begin
    A = Alphabet([:x, :X], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    RF = LowCohomologySOS.suitable_group_ring([a^2])
    b = RF.basis
    mstr = RF.mstructure
    @test RF.object == F
    @test Set(collect(b)) ==
          Set([one(F), a, a^(-1), a^2, a^(-2), a^3, a^(-3), a^4, a^(-4)])
    @test mstr[getindex(b, a), getindex(b, a)] == getindex(b, a^2)
    @test mstr[getindex(b, a), getindex(b, a^2)] == getindex(b, a^3)
    @test mstr[getindex(b, a^2), getindex(b, a^(-1))] == getindex(b, a)

    A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    F2 = FreeGroup(A)
    x, y = Groups.gens(F2)
    RF2 = LowCohomologySOS.suitable_group_ring([x, y])
    b2 = RF2.basis
    mstr2 = RF2.mstructure
    @test RF2.object == F2
    @test Set(collect(b2)) == Set([
        one(F2),
        x,
        x^(-1),
        y,
        y^(-1),
        x^2,
        x^(-2),
        y^2,
        y^(-2),
        x * y,
        x * y^(-1),
        x^(-1) * y,
        x^(-1) * y^(-1),
        y * x,
        y * x^(-1),
        y^(-1) * x,
        y^(-1) * x^(-1),
    ])
    @test mstr2[getindex(b2, x), getindex(b2, x)] == getindex(b2, x^2)
    @test mstr2[getindex(b2, y), getindex(b2, y^(-1))] == getindex(b2, one(F2))
    @test mstr2[getindex(b2, x^(-1)), getindex(b2, y)] ==
          getindex(b2, x^(-1) * y)
end
