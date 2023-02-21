@testset "sln embeddings" begin
    i = LowCohomologySOS.sln_slm_embedding(3,4)

    sl3 = i.source
    sl4 = i.target

    e12, e13, e21, e23, e31, e32 = gens(sl3)
    E12, E13, E14, E21, E23, E24, E31, E32, E34, E41, E42, E43 = gens(sl4)

    @test i(e12) == E12
    @test i(e13) == E13
    @test i(e21) == E21
    @test i(e23) == E23
    @test i(e31) == E31
    @test i(e32) == E32

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