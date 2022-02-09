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

    @testset "one homomorphism" begin
        one_hom = PropertyT_new.Homomorphism(
            (i,source, target) -> word(one(target)),
            F,
            G
        )
        test_homomorphism(one_hom)
        @test isone(one_hom(x))
        @test isone(one_hom(y))

        @test LowCohomologySOS.embed(one_hom, zero(RF), RG) == zero(RG)
        @test LowCohomologySOS.embed(one_hom, RF(x), RG) == one(RG)
        @test LowCohomologySOS.embed(one_hom, RF(y), RG) == one(RG)
        @test LowCohomologySOS.embed(one_hom, RF(x * y) + RF(y^2), RG) == 2 * one(RG)
    end

    @testset "quotient homomorphism" begin
        quotient_hom = let source = F, target = G
            function f(i, source, target)
                if alphabet(source) == alphabet(target)
                    Groups.word_type(target)([i])
                else
                    throw("Unsupported")
                end
            end
            PropertyT_new.Homomorphism(f, source, target)
        end

        test_homomorphism(quotient_hom)

        @test !isone(quotient_hom(x))
        @test !isone(quotient_hom(y))
        @test quotient_hom(x*y) == quotient_hom(y*x)
        @test isone(quotient_hom(Groups.commutator(x,y)))

        @test LowCohomologySOS.embed(quotient_hom, zero(RF), RG) == zero(RG)
        let (x, y) = Groups.gens(F)
            xx = quotient_hom(x)
            yy = quotient_hom(y)
            @test LowCohomologySOS.embed(quotient_hom, RF(x), RG) == RG(xx)
            @test LowCohomologySOS.embed(quotient_hom, RF(y), RG) == RG(yy)
            @test LowCohomologySOS.embed(quotient_hom, RF(x * y) + RF(y^2), RG) == RG(xx * yy) + RG(yy^2)
        end
    end
end

@testset "Fox derivatives" begin
    @testset "via group ring of the free group" begin
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

    @testset "via FreeGroup, return vectors" begin
        F = FreeGroup(4)

        c, g = fox_derivative(F, one(F), 1)
        @test isempty(c) && isempty(g)
        c, g  = fox_derivative(F, gens(F, 1), 1)
        @test c == [1] && g == [one(F)]

        c, g = fox_derivative(F, inv(gens(F, 1)), 1)
        @test c == [-1] && g == [inv(gens(F, 1))]
        c, g = fox_derivative(F, inv(gens(F, 1)), 2)
        @test isempty(c) && isempty(g)

        a,b,c,d = gens(F)
        u = a*b^-1*c*d^-1

        cf, elts = fox_derivative(F, u, 1)
        @test cf == [1] && elts == [one(F)]

        cf, elts = fox_derivative(F, u, 2)
        @test cf == [-1] && elts == [a*b^-1]

        cf, elts = fox_derivative(F, u, 3)
        @test cf == [1] && elts == [a*b^-1]

        cf, elts = fox_derivative(F, u, 4)
        @test cf == [-1] && elts == [a*b^-1*c*d^-1]
    end

    @testset "agreement of two definitions" begin
        F = FreeGroup(4)
        RF = LowCohomologySOS.group_ring(F, 3)

        A = alphabet(F)
        for u in [F(rand(1:length(A), 6)) for _ in 1:4]
            for i in 1:4
                (coeffs, elts) = fox_derivative(F, u, i)

                fd = if isempty(coeffs)
                    @test isempty(elts)
                    zero(RF)
                else
                    sum(c*RF(g) for (c, g) in zip(coeffs, elts))
                end

                try
                    @test LowCohomologySOS.fox_derivative(RF, u, i) == fd
                catch err
                    @error "failed (?) test for" u, i
                    rethrow(err)
                end
            end
        end
    end
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

# TODO: tests for suitable_group_ring
