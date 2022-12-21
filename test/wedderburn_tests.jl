@testset "conjugation action on Nielsen generators of SAut(Fₙ)" begin
    _conj = LowCohomologySOS._conj

    S₃ = PermutationGroups.SymmetricGroup(3)
    for σ ∈ S₃
        i, j = rand(1:3), rand(1:3)
        @test _conj(Groups.Transvection(:ϱ, i, j), σ) == Groups.Transvection(:ϱ, σ[i], σ[j])
        @test _conj(Groups.Transvection(:λ, i, j), σ) == Groups.Transvection(:λ, σ[i], σ[j])
        @test _conj(Groups.Transvection(:ϱ, i, j, true), σ) == Groups.Transvection(:ϱ, σ[i], σ[j], true)
        @test _conj(Groups.Transvection(:λ, i, j, true), σ) == Groups.Transvection(:λ, σ[i], σ[j], true)
    end

    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(3))
    for g in Σ
        i, j = rand(1:3), rand(1:3)
        if isone(g.n.elts[i])
            @test _conj(Groups.Transvection(:ϱ, i, j), g) == Groups.Transvection(:ϱ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? false : true)
            @test _conj(Groups.Transvection(:λ, i, j), g) == Groups.Transvection(:λ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? false : true)
            @test _conj(Groups.Transvection(:ϱ, i, j, true), g) == Groups.Transvection(:ϱ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? true : false)
            @test _conj(Groups.Transvection(:λ, i, j, true), g) == Groups.Transvection(:λ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? true : false)
        else
            @test _conj(Groups.Transvection(:ϱ, i, j), g) == Groups.Transvection(:λ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? true : false)
            @test _conj(Groups.Transvection(:λ, i, j), g) == Groups.Transvection(:ϱ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? true : false)
            @test _conj(Groups.Transvection(:ϱ, i, j, true), g) == Groups.Transvection(:λ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? false : true)
            @test _conj(Groups.Transvection(:λ, i, j, true), g) == Groups.Transvection(:ϱ, inv(g.p)[i], inv(g.p)[j], isone(g.n.elts[j]) ? false : true)
        end
    end
end

@testset "permutations of Nielsen transvections corresponding to wreath product elements" begin
    SAutF₃ = Groups.SpecialAutomorphismGroup(FreeGroup(3))
    SAutF₃_alphabet = alphabet(SAutF₃)
    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(3))
    _conj = LowCohomologySOS._conj
    action = LowCohomologySOS.AlphabetPermutation(SAutF₃_alphabet, Σ, _conj)
    for σ ∈ Σ
        σ_perm = action.perms[σ]
        for l in SAutF₃_alphabet.letters
            @test σ_perm[SAutF₃_alphabet[l]] == SAutF₃_alphabet[_conj(l, σ)]
        end
    end
end

@testset "constructions of basis elements for the matrix SOS problem" begin
    SAutF₂ = Groups.SpecialAutomorphismGroup(FreeGroup(2))
    S = let s = Groups.gens(SAutF₂)
        [s; inv.(s)]
    end
    S = unique!(S)
    half_basis = S
    basis, sizes = Groups.wlmetric_ball(S, radius = 2)

    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)

    @test length(constraints_basis) == length(S)^2*length(basis)
    @test length(psd_basis) == length(S)*length(half_basis)

    gl = length(S)
    hbl = sizes[1]
    for l in 1:3
        i, k = rand(1:gl), rand(1:hbl)
        psd_basis_element = LowCohomologySOS.PSDBasisElement(S[i], basis[k])
        @test psd_basis_element isa LowCohomologySOS.PSDBasisElement
        @test psd_basis_element.generator isa Groups.GroupElement
        @test psd_basis_element.basis_elt isa Groups.GroupElement
        @test psd_basis_element.generator == S[i]
        @test psd_basis_element.basis_elt == basis[k]
        @test psd_basis_element == LowCohomologySOS.PSDBasisElement(S[i], basis[k])
    end
end

@testset "action of wreath product on basis elements" begin
    SAutF₃ = Groups.SpecialAutomorphismGroup(FreeGroup(3))
    S = let s = Groups.gens(SAutF₃)
        [s; inv.(s)]
    end
    S = unique!(S)
    gl = length(S)
    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(3))
    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)

    for l in 1:3
        i, k = rand(1:length(S)), rand(1:length(basis))
        basis_element = LowCohomologySOS.PSDBasisElement(S[i], basis[k])
        for σ ∈ Σ
            be_after_action = SymbolicWedderburn.action(actions, σ, basis_element)
            @test be_after_action == LowCohomologySOS.PSDBasisElement(
                SAutF₃(word(basis_element.generator)^actions.alphabet_perm.perms[σ]),
                SAutF₃(word(basis_element.basis_elt)^actions.alphabet_perm.perms[σ])
            )
        end
    end

    C_ind = CartesianIndices((length(basis),length(S),length(S)))
    L_ind = LinearIndices(C_ind)
    for l in 1:3
        basis_element = LowCohomologySOS.LinIdx(rand(1:length(S)^2*length(basis)))
        be_k, be_j, be_i = Tuple(C_ind[basis_element.id])
        for σ ∈ Σ
            be_after_action = SymbolicWedderburn.action(actions, σ, basis_element)
            be_k_σ = be_k^SymbolicWedderburn.induce(actions.basis_action, σ)
            be_j_σ = be_j^SymbolicWedderburn.induce(actions.S_action, σ)
            be_i_σ = be_i^SymbolicWedderburn.induce(actions.S_action, σ)
            rhs = LowCohomologySOS.LinIdx(L_ind[be_k_σ,be_j_σ,be_i_σ])
            @test be_after_action == rhs
        end
    end
end

@testset "wedderburn_decomposition_matrix" begin
    SAutF₂ = Groups.SpecialAutomorphismGroup(FreeGroup(2))
    S = let s = Groups.gens(SAutF₂)
        [s; inv.(s)]
    end
    S = unique!(S)
    half_basis = S
    basis, sizes = Groups.wlmetric_ball(S, radius = 2)
    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(2))
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)

    @test eltype(w_dec_matrix.basis) <: LowCohomologySOS.LinIdx
    @test length(w_dec_matrix.basis) == length(S)^2*length(basis)
    @test length(w_dec_matrix.invariants[rand(1:456)]) == length(S)^2*length(basis)
    for i in 1:length(w_dec_matrix.Uπs)
        @test size(w_dec_matrix.Uπs[i])[2] == length(S)*length(half_basis)
    end

    sa = SymbolicWedderburn.symmetry_adapted_basis(Rational{Int}, Σ, actions, psd_basis)
    mπs = SymbolicWedderburn.multiplicity.(sa)

    for i in 1:length(mπs)
        @test size(w_dec_matrix.Uπs[i])[1] == mπs[i]
    end
end

@testset "expressing bases: idies and triples agreement and proper order" begin
    bs, n = Int32(rand(10:1_700_000)), Int32(rand(2:49))

    C_indicies = CartesianIndices((bs, n, n))
    L_indices = LinearIndices(C_indicies)

    for it in 1:100
        k, j, i = rand(UInt32(1):bs), rand(UInt32(1):n), rand(UInt32(1):n)
        @test Tuple(C_indicies[L_indices[k,j,i]]) == (k,j,i)
    end
    for it in 1:100
        id = rand(UInt32(1):n^2*bs)
        k,j,i = Tuple(C_indicies[id])
        @test L_indices[k,j,i] == id
    end
end