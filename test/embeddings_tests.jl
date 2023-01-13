SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)

@testset "sln embeddings" begin
    sl3, sl4 = SL(3,Int8), SL(4,Int8)

    Σ = PermutationGroups.SymmetricGroup(4)

    e12, e13, e21, e23, e31, e32 = gens(sl3)
    _idx(k) = ((i,j) for i in 1:k for j in 1:k if i≠j)
    E(i,j) = MatrixGroups.ElementaryMatrix{4}(i,j, Int8(1))
    gens_from_elmatrices = Dict(
         E(i,j) => sl4([alphabet(sl4)[E(i,j)]]) for (i,j) in _idx(4)
    )
    for σ in Σ
        i = LowCohomologySOS.sln_slm_embedding(sl3, sl4, σ)
        @test i(e12) == gens_from_elmatrices[E(1^σ,2^σ)]
        @test i(e13) == gens_from_elmatrices[E(1^σ,3^σ)]
        @test i(e21) == gens_from_elmatrices[E(2^σ,1^σ)]
        @test i(e23) == gens_from_elmatrices[E(2^σ,3^σ)]
        @test i(e31) == gens_from_elmatrices[E(3^σ,1^σ)]
        @test i(e32) == gens_from_elmatrices[E(3^σ,2^σ)]
    end

    idx = one(Σ)

    i = LowCohomologySOS.sln_slm_embedding(sl3, sl4, idx)

    E12, E13, E14, E21, E23, E24, E31, E32, E34, E41, E42, E43 = gens(sl4)

    RG  = LowCohomologySOS.group_ring(sl3, 2)
    M = [
        one(RG) one(RG) zero(RG) zero(RG) RG(e12) zero(RG);
        RG(e13) one(RG) one(RG) zero(RG) zero(RG) zero(RG);
        one(RG) one(RG) RG(e21) zero(RG) zero(RG) zero(RG);
        one(RG) one(RG) zero(RG) zero(RG) RG(e23) zero(RG);
        one(RG) one(RG) zero(RG) zero(RG) RG(e31) zero(RG);
        one(RG) one(RG) zero(RG) zero(RG) RG(e32) zero(RG)
    ]

    RG_prime = LowCohomologySOS.group_ring(sl4,2)
    M_emb = LowCohomologySOS.embed_matrix(M, i, RG_prime)
    _0 = zero(RG_prime)
    M_emb_proper = [
        one(RG_prime) one(RG_prime) _0 _0 _0 _0 RG_prime(E12) _0 _0 _0 _0 _0;
        RG_prime(E13) one(RG_prime) _0 one(RG_prime) _0 _0 _0 _0 _0 _0 _0 _0;
        _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0;
        one(RG_prime) one(RG_prime) _0 RG_prime(E21) _0 _0 _0 _0 _0 _0 _0 _0;
        one(RG_prime) one(RG_prime) _0 _0 _0 _0 RG_prime(E23) _0 _0 _0 _0 _0;
        _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0;
        one(RG_prime) one(RG_prime) _0 _0 _0 _0 RG_prime(E31) _0 _0 _0 _0 _0;
        one(RG_prime) one(RG_prime) _0 _0 _0 _0 RG_prime(E32) _0 _0 _0 _0 _0;
        _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0;
        _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0;
        _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0;
        _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0 _0
    ]
    
    @test M_emb == M_emb_proper
end