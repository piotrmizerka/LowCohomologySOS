@testset "SL(3, ℤ)" begin
    A = Alphabet(
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
    F_sl_3_z = FreeGroup(A)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F_sl_3_z)

    N = 3
    SL(n, R) = PropertyT_new.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    h = let source = F_sl_3_z, target = SL₃ℤ
        function hom(i, source, target::PropertyT_new.SpecialLinearGroup{N, T}) where {N, T}
            As = alphabet(source)
            s = String(As[i])
            (i,j,isinv) = let m = match(r"^(E|e)(\d)(\d)$", s)
                @assert !isnothing(m)
                isinv = isuppercase(Char(first(m[1])))
                i = parse(Int, m[2])
                j = parse(Int, m[3])
                i,j, isinv
            end

            At = alphabet(target)
            eij = PropertyT_new.MatrixGroups.ElementaryMatrix{N}(i, j, T(isinv ? -1 : 1))
            return Groups.word_type(target)([At[eij]])
        end
        PropertyT_new.Homomorphism(hom, source, target)
    end

    @testset "homomorphism basic tests" begin
        test_homomorphism(h)
    end
    
    @testset "check if relations are well-defined and preserved by the quotient homomorphism" begin
        @test h(Groups.commutator(e12, e13)) == one(SL₃ℤ)
        @test h(Groups.commutator(e12, e32)) == one(SL₃ℤ)
        @test h(Groups.commutator(e13, e23)) == one(SL₃ℤ)
        @test h(Groups.commutator(e23, e21)) == one(SL₃ℤ)
        @test h(Groups.commutator(e21, e31)) == one(SL₃ℤ)
        @test h(Groups.commutator(e31, e32)) == one(SL₃ℤ)
        @test h(Groups.commutator(e12, e23)) == h(e13)
        @test h(Groups.commutator(e13, e32)) == h(e12)
        @test h(Groups.commutator(e21, e13)) == h(e23)
        @test h(Groups.commutator(e23, e31)) == h(e21)
        @test h(Groups.commutator(e31, e12)) == h(e32)
        @test h(Groups.commutator(e32, e21)) == h(e31)
        @test h((e12*e21^(-1)*e12)^4) == one(SL₃ℤ) # Gersten relation
    end
    

    _idx(n) = ((i,j) for i in 1:n for j in 1:n if i≠j)

    S = Dict((i,j) =>
        let eij = PropertyT_new.MatrixGroups.ElementaryMatrix{N}(i,j, Int8(1))
            SL₃ℤ([alphabet(SL₃ℤ)[eij]])
        end
        for (i,j) in _idx(N)
    )

    S_free = Dict(
        (i,j) => (s = Symbol('e', i, j); F = F_sl_3_z; F([alphabet(F)[s]]))
        for (i,j) in _idx(N)
    )

    @testset "homomorphism specific tests" begin
        @test all(h(S_free[(i,j)]) == S[(i,j)] for (i,j) in _idx(N))
        for (i,j) in _idx(N)
            for (k,l) in _idx(N)
                g_free = S_free[(i,j)] * S_free[(k,l)]
                g = S[(i,j)] * S[(k,l)]
                @test h(g_free) == g
            end
        end
    end

    @testset "group ring embedding" begin
        RF_sl_3_z = LowCohomologySOS.group_ring(F_sl_3_z, 1)

        RSL₃ℤ, sizes = let G = SL₃ℤ, S = union!(collect(values(S)), inv.(values(S)))
            half_basis, sizes = Groups.wlmetric_ball(S, radius = 1)
            LowCohomologySOS.group_ring(G, half_basis, true), sizes
        end

        @test LowCohomologySOS.embed(h, RF_sl_3_z(S_free[(1,2)] * S_free[(1,3)]^-1), RSL₃ℤ) ==
              RSL₃ℤ(S[(1,2)] * S[(1,3)]^(-1))

        @test all(_idx(N)) do (i,j)
            g_free, g = S_free[(i,j)], S[(i,j)]
            LowCohomologySOS.embed(h, RF_sl_3_z(g_free), RSL₃ℤ) == RSL₃ℤ(g)
        end
    end
end