@testset "invariant_constraint_matrix" begin
    SAutF₂ = Groups.SpecialAutomorphismGroup(FreeGroup(2))
    S = let s = Groups.gens(SAutF₂)
        [s; inv.(s)]
    end
    S = unique!(S)
    half_basis = S
    basis, sizes = Groups.wlmetric_ball(S, radius = 2)
    _conj = LowCohomologySOS._conj
    ℝSAutF₂_star = LowCohomologySOS.group_ring(SAutF₂, half_basis, star_multiplication = true)
    cnstrs = LowCohomologySOS.constraints(ℝSAutF₂_star.mstructure)

    i, j, e = rand(1:length(S)), rand(1:length(S)), rand(1:length(basis))

    orbit_representative = LowCohomologySOS.TensorSupportElement(S[i], S[j], basis[e])
    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(2))
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, _conj)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    A_gs = Dict(g => A_g for (A_g, g) in zip(cnstrs, basis))

    inv_cnstr_matrix = LowCohomologySOS.invariant_constraint_matrix(
        orbit_representative,
        Σ,
        action,
        psd_basis,
        A_gs,
        length(S)
    )

    @test size(inv_cnstr_matrix) == (length(psd_basis), length(psd_basis))
    @info inv_cnstr_matrix
    # check invariance of inv_cnstr_matrix: this test is failed - study the invariance!
    for i in 1:length(S)
        for j in 1:length(S)
            inv_cnstr_matrix_i_j = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(i, j)]
            for σ ∈ Σ
                σ_tse = SymbolicWedderburn.action(action, σ, LowCohomologySOS.TensorSupportElement(S[i], S[j], first(S)))
                σi, σj = word(σ_tse.row_generator)[1], word(σ_tse.column_generator)[1]
                inv_cnstr_matrix_σi_σj = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(σi, σj)]
                for l in 1:3
                    e, f = rand(1:length(half_basis)), rand(1:length(half_basis))
                    σe, σf = (word(half_basis[e])^action.perms[σ])[1], (word(half_basis[f])^action.perms[σ])[1]
                    @test inv_cnstr_matrix_σi_σj[σe,σf] == inv_cnstr_matrix_i_j[e,f]
                end
            end
        end
    end

    M = [i ≠ j ? zero(ℝSAutF₂_star) : one(ℝSAutF₂_star)+ℝSAutF₂_star(S[2]) for i in 1:length(S), j in 1:length(S)]
    order_unit = [i ≠ j ? zero(ℝSAutF₂_star) : one(ℝSAutF₂_star) for i in 1:length(S), j in 1:length(S)]
    w_dec_matrix = LowCohomologySOS.wedderburn_decomposition_matrix(Σ, basis, half_basis, S)

    sos_pr_sym = LowCohomologySOS.sos_problem_symmetrized(
        M,
        order_unit,
        w_dec_matrix,
        Σ,
        psd_basis,
        S
    )

    # TODO: more tests - especially the missing ones for sos_problem_symmetrized
end