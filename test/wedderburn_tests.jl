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
        if isone(g.n.elts[g.p[i]])
            @test _conj(Groups.Transvection(:ϱ, i, j), g) == Groups.Transvection(:ϱ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? false : true)
            @test _conj(Groups.Transvection(:λ, i, j), g) == Groups.Transvection(:λ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? false : true)
            @test _conj(Groups.Transvection(:ϱ, i, j, true), g) == Groups.Transvection(:ϱ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? true : false)
            @test _conj(Groups.Transvection(:λ, i, j, true), g) == Groups.Transvection(:λ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? true : false)
        else
            @test _conj(Groups.Transvection(:ϱ, i, j), g) == Groups.Transvection(:λ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? true : false)
            @test _conj(Groups.Transvection(:λ, i, j), g) == Groups.Transvection(:ϱ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? true : false)
            @test _conj(Groups.Transvection(:ϱ, i, j, true), g) == Groups.Transvection(:λ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? false : true)
            @test _conj(Groups.Transvection(:λ, i, j, true), g) == Groups.Transvection(:ϱ, g.p[i], g.p[j], isone(g.n.elts[g.p[j]]) ? false : true)
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
    bl = length(basis)
    for i in 1:3
        row_gen, col_gen, entry = rand(1:gl), rand(1:gl), rand(1:bl)
        basis_element = LowCohomologySOS.TensorSupportElement(S[row_gen], S[col_gen], basis[entry])
        @test basis_element isa LowCohomologySOS.TensorSupportElement
        @test basis_element.row_generator isa Groups.GroupElement
        @test basis_element.column_generator isa Groups.GroupElement
        @test basis_element.entry isa Groups.GroupElement
        @test basis_element.row_generator == S[row_gen]
        @test basis_element.column_generator == S[col_gen]
        @test basis_element.entry == basis[entry]
        @test basis_element == LowCohomologySOS.TensorSupportElement(S[row_gen], S[col_gen], basis[entry])
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
    _conj = LowCohomologySOS._conj
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, _conj)

    for i in 1:3
        row_gen, col_gen, entry = rand(1:gl), rand(1:gl), rand(1:gl)
        basis_element = LowCohomologySOS.TensorSupportElement(S[row_gen], S[col_gen], S[entry])
        for σ ∈ Σ
            be_after_action = SymbolicWedderburn.action(action, σ, basis_element)
            @test be_after_action == LowCohomologySOS.TensorSupportElement(
                SAutF₃(word(basis_element.row_generator)^action.perms[σ]),
                SAutF₃(word(basis_element.column_generator)^action.perms[σ]),
                SAutF₃(word(basis_element.entry)^action.perms[σ])
            )
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
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)

    w_dec_matrix = LowCohomologySOS.wedderburn_decomposition_matrix(Σ, constraints_basis, psd_basis, S)

    @test eltype(w_dec_matrix.basis) == LowCohomologySOS.TensorSupportElement
    @test length(w_dec_matrix.basis) == length(S)^2*length(basis)
    @test length(w_dec_matrix.invariants[rand(1:456)]) == length(S)^2*length(basis)
    for i in 1:length(w_dec_matrix.Uπs)
        @test size(w_dec_matrix.Uπs[i])[2] == length(S)*length(half_basis)
    end

    _conj = LowCohomologySOS._conj
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, _conj)
    sa = SymbolicWedderburn.symmetry_adapted_basis(Rational{Int}, Σ, action, psd_basis)
    mπs = SymbolicWedderburn.multiplicity.(sa)

    for i in 1:length(mπs)
        @test size(w_dec_matrix.Uπs[i])[1] == mπs[i]
    end
end