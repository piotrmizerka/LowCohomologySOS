SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)
SAut_F(n) = Groups.SpecialAutomorphismGroup(FreeGroup(n))

@testset "relations for SL(n,ℤ) and SAut(Fₙ)" begin

    @testset "symmetric group actions" begin
        symmetric_action = true

        @testset "SL(n,ℤ)" begin
            G = SL(3,Int8)
            S = gens(G)
            F_G = Groups.FreeGroup(alphabet(G))

            relations = LowCohomologySOS.relations(G, F_G, S, symmetric_action, 3)

            e12, e13, e21, e23, e31, e32 = gens(F_G)

            @test Set(relations) == Set([
                e12 * e13 * e12^(-1) * e13^(-1),
                e13 * e12 * e13^(-1) * e12^(-1),
                e21 * e23 * e21^(-1) * e23^(-1),
                e23 * e21 * e23^(-1) * e21^(-1),
                e31 * e32 * e31^(-1) * e32^(-1),
                e32 * e31 * e32^(-1) * e31^(-1),
                e21 * e31 * e21^(-1) * e31^(-1),
                e31 * e21 * e31^(-1) * e21^(-1),
                e12 * e32 * e12^(-1) * e32^(-1),
                e32 * e12 * e32^(-1) * e12^(-1),
                e13 * e23 * e13^(-1) * e23^(-1),
                e23 * e13 * e23^(-1) * e13^(-1),
                e12 * e23 * e12^(-1) * e23^(-1) * e13^(-1),
                e13 * e32 * e13^(-1) * e32^(-1) * e12^(-1),
                e21 * e13 * e21^(-1) * e13^(-1) * e23^(-1),
                e23 * e31 * e23^(-1) * e31^(-1) * e21^(-1),
                e31 * e12 * e31^(-1) * e12^(-1) * e32^(-1),
                e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1)
            ])
        end

        @testset "SAut(Fₙ)" begin
            G = SAut_F(3)
            S = gens(G)
            F_G = Groups.FreeGroup(alphabet(G))

            relations = LowCohomologySOS.relations(G, F_G, S, symmetric_action, 3)

            r12, r13, r21, r23, r31, r32, l12, l13, l21, l23, l31, l32 = gens(F_G)

            quadruple_relations = [
                r32^(-1) * r12^(-1) * r32 * r12,
                r23^(-1) * r13^(-1) * r23 * r13,
                r31^(-1) * r21^(-1) * r31 * r21,
                r13^(-1) * r23^(-1) * r13 * r23,
                r21^(-1) * r31^(-1) * r21 * r31,
                r12^(-1) * r32^(-1) * r12 * r32,
                l32^(-1) * l12^(-1) * l32 * l12,
                l23^(-1) * l13^(-1) * l23 * l13,
                l31^(-1) * l21^(-1) * l31 * l21,
                l13^(-1) * l23^(-1) * l13 * l23,
                l21^(-1) * l31^(-1) * l21 * l31,
                l12^(-1) * l32^(-1) * l12 * l32,
                l32^(-1) * r12^(-1) * l32 * r12,
                l23^(-1) * r13^(-1) * l23 * r13,
                l31^(-1) * r21^(-1) * l31 * r21,
                l13^(-1) * r23^(-1) * l13 * r23,
                l21^(-1) * r31^(-1) * l21 * r31,
                l12^(-1) * r32^(-1) * l12 * r32,
                l12^(-1) * r12^(-1) * l12 * r12,
                l13^(-1) * r13^(-1) * l13 * r13,
                l21^(-1) * r21^(-1) * l21 * r21,
                l23^(-1) * r23^(-1) * l23 * r23,
                l31^(-1) * r31^(-1) * l31 * r31,
                l31^(-1) * r31^(-1) * l31 * r31
            ]

            @test issubset(quadruple_relations, relations)
            @test length(relations) == 78
        end
    end

    @testset "wreath product actions" begin
        symmetric_action = false

        @testset "SL(n,ℤ)" begin
            G = SL(3,Int8)
            S = let s = gens(G)
                [s; inv.(s)]
            end
            F_G = Groups.FreeGroup(length(S))

            relations = LowCohomologySOS.relations(G, F_G, S, symmetric_action, 3)

            @test length(relations) == 4*18+6*2+6*4
        end

        @testset "SAut(Fₙ)" begin
            G = SAut_F(3)
            S = let s = gens(G)
                [s; inv.(s)]
            end
            F_G = Groups.FreeGroup(length(S))

            relations = LowCohomologySOS.relations(G, F_G, S, symmetric_action, 3)

            pairs_rels_number = 24
            quadruples_1_rels_number = 2*6*4
            quadruples_2_rels_number = 2*(6+6+6)*4
            triples_rels_number = 8*6
            @test length(relations) == pairs_rels_number+
                                        quadruples_1_rels_number+
                                        quadruples_2_rels_number+
                                        triples_rels_number
        end
    end
end

@testset "star conjugation kills yields the squared Laplacian" begin

    function star_conj_test(G, wreath)
        S = gens(G)
        S_inv = [S; inv.(S)]
        half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2)
        proper_S = wreath ? S_inv : S

        Δ₁, Iₙ, Δ₁⁺, Δ₁⁻, d₀_ = LowCohomologySOS.laplacians(G, half_basis, proper_S, twist_coeffs = false);

        @test d₀_'*Δ₁*d₀_ == d₀_'*d₀_*d₀_'*d₀_
    end

    @testset "symmetric group actions" begin

        @testset "SL(n,ℤ)" begin
            star_conj_test(SL(3,Int8), false)
        end

        @testset "SAut(Fₙ)" begin
            star_conj_test(SAut_F(3), false)
        end

    end

    @testset "wreath product actions" begin

        @testset "SL(n,ℤ)" begin
            star_conj_test(SL(3,Int8), true)
        end

        @testset "SAut(Fₙ)" begin
            star_conj_test(SAut_F(3), true)
        end

    end
end
