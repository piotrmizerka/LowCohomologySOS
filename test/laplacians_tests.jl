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

            commutator_relations = [
                r32 * r12 * r32^(-1) * r12^(-1),
                r23 * r13 * r23^(-1) * r13^(-1),
                r31 * r21 * r31^(-1) * r21^(-1),
                r13 * r23 * r13^(-1) * r23^(-1),
                r21 * r31 * r21^(-1) * r31^(-1),
                r12 * r32 * r12^(-1) * r32^(-1),
                l32 * l12 * l32^(-1) * l12^(-1),
                l23 * l13 * l23^(-1) * l13^(-1),
                l31 * l21 * l31^(-1) * l21^(-1),
                l13 * l23 * l13^(-1) * l23^(-1),
                l21 * l31 * l21^(-1) * l31^(-1),
                l12 * l32 * l12^(-1) * l32^(-1),
                l32 * r12 * l32^(-1) * r12^(-1),
                l23 * r13 * l23^(-1) * r13^(-1),
                l31 * r21 * l31^(-1) * r21^(-1),
                l13 * r23 * l13^(-1) * r23^(-1),
                l21 * r31 * l21^(-1) * r31^(-1),
                l12 * r32 * l12^(-1) * r32^(-1),
                l12 * r12 * l12^(-1) * r12^(-1),
                l13 * r13 * l13^(-1) * r13^(-1),
                l21 * r21 * l21^(-1) * r21^(-1),
                l23 * r23 * l23^(-1) * r23^(-1),
                l31 * r31 * l31^(-1) * r31^(-1),
                l31 * r31 * l31^(-1) * r31^(-1)
            ]

            @test issubset(commutator_relations, relations)
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

@testset "star conjugation yields the squared Laplacian" begin

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

@testset "Adj⁺ for SLₙ(ℤ)" begin
    N = 3
    G = SL(N,Int8)
    A = alphabet(G)
    eij_id = Dict(A[i] => i for i in eachindex(A))
    e(i,j) = G([eij_id[MatrixGroups.ElementaryMatrix{N}(i,j,Int8(1))]])
    S = gens(G)
    S_inv = [S;inv.(S)]
    half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2)
    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(
        G, half_basis, S, sq_adj_op_ = "adj", twist_coeffs = false)
    RG = parent(first(Δ₁⁺))

    function ijij(i,j)
        range_as_list = [i for i in 1:N]
        range_no_ij = deleteat!(copy(range_as_list), findall(l->l∈[i,j],copy(range_as_list)))

        result = 2*sum(
            (one(RG)-RG(e(i,l)))'*(one(RG)-RG(e(i,l)))+
            (one(RG)-RG(e(l,j)))'*(one(RG)-RG(e(l,j)))
        for l in range_no_ij)+
        sum(
            (one(RG)-RG(e(i,l)*e(j,l)))'*(one(RG)-RG(e(i,l)*e(j,l)))+
            (RG(e(l,i))-RG(e(l,j)))'*(RG(e(l,i))-RG(e(l,j)))
        for l in range_no_ij)+
        (N-2)*one(RG)

        return result
    end
    function ijik(i,j,k)
        return 2*(one(RG)-RG(e(i,k)))'*(RG(e(i,j))-one(RG))-(one(RG)-RG(e(i,k)*e(j,k)))'-(one(RG)-RG(e(i,j)*e(k,j)))
    end
    function ijki(i,j,k)
        return (RG(e(k,i))-RG(e(k,j)))'*(one(RG)-RG(e(k,j)*e(i,j)))
    end
    function ijjk(i,j,k)
        return (one(RG)-RG(e(i,k)*e(j,k)))'*(RG(e(i,j))-RG(e(i,k)))
    end
    function ijkj(i,j,k)
        return 2*(one(RG)-RG(e(k,j)))'*(RG(e(i,j))-one(RG))-(RG(e(k,i))-RG(e(k,j)))'-(RG(e(i,k))-RG(e(i,j)))
    end

    gen_ij = Dict(s => (A[word(s)[1]].i,A[word(s)[1]].j) for s in S)
    for s in eachindex(S)
        i,j = gen_ij[S[s]]
        for t in eachindex(S)
            k,l = gen_ij[S[t]]
            if i == k && j == l
                @test Δ₁⁺[s,t] == ijij(i,j)
            elseif length(collect(Set([i,j,k,l]))) == 3
                if i == k
                    @test Δ₁⁺[s,t] == ijik(i,j,l)
                elseif i == l
                    @test Δ₁⁺[s,t] == ijki(i,j,k)
                elseif j == k
                    @test Δ₁⁺[s,t] == ijjk(i,j,l)
                elseif j == l
                    @test Δ₁⁺[s,t] == ijkj(i,j,k)
                end
            else
                @test Δ₁⁺[s,t] == zero(RG)
            end
        end
    end
end