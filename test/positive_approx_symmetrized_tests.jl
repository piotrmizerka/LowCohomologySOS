using SparseArrays
using Kronecker

function act_on_matrix(
    P, 
    σ::Groups.GroupElement,
    half_basis,
    half_basis_idies::Dict,
    action
)
    result = zeros(eltype(P), size(P)[1], size(P)[2])

    for e in half_basis
        for f in half_basis
            id_e, id_f = half_basis_idies[e], half_basis_idies[f]
            eσ_inv, fσ_inv = SymbolicWedderburn.action(action, σ^(-1), e), SymbolicWedderburn.action(action, σ^(-1), f)
            id_eσ_inv, id_fσ_inv = half_basis_idies[eσ_inv], half_basis_idies[fσ_inv]
            result[id_e, id_f] = P[id_eσ_inv, id_fσ_inv]
        end
    end

    return result
end

# @testset "symmetrized matrix SOS problem" begin
    N = 4
    half_radius = 1

    SAutF_N = Groups.SpecialAutomorphismGroup(FreeGroup(N))
    S = let s = Groups.gens(SAutF_N)
        [s; inv.(s)]
    end
    S = unique!(S)
    basis, sizes = Groups.wlmetric_ball(S, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]
    ℝSAutF_N_star = LowCohomologySOS.group_ring(SAutF_N, half_basis, star_multiplication = true)
    cnstrs =  dropzeros!.(sparse.(LowCohomologySOS.constraints(ℝSAutF_N_star.mstructure)))
    A_gs_cart = [findall(!iszero,c) for c in cnstrs]

    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(N))
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(ℝSAutF_N_star.basis, half_basis, S)
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)

    inv_cnstr_matrix = LowCohomologySOS.invariant_constraint_matrix(
        rand(w_dec_matrix.invariants),
        A_gs_cart,
        length(half_basis),
        length(collect(Σ))
    )

    @test size(inv_cnstr_matrix) == (length(psd_basis), length(psd_basis))

    half_basis_idies = Dict(half_basis[i] => i for i in 1:length(half_basis))
    S_idies = Dict(S[i] => i for i in 1:length(S))
    for i in 1:length(S)
        for j in 1:length(S)
            lhs = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(i, j)]
            for σ ∈ Σ
                tse_σ_inv = SymbolicWedderburn.action(action, σ^(-1), LowCohomologySOS.TensorSupportElement(S[i], S[j], first(S)))
                iσ_inv, jσ_inv = S_idies[tse_σ_inv.i], S_idies[tse_σ_inv.j]
                P = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(iσ_inv, jσ_inv)]
                rhs = act_on_matrix(P, σ, half_basis, half_basis_idies, action)
                @test lhs == rhs
            end
        end
    end

    M = [i ≠ j ? zero(ℝSAutF_N_star) : one(ℝSAutF_N_star)+ℝSAutF_N_star(S[2]) for i in 1:length(S), j in 1:length(S)]
    order_unit = [i ≠ j ? zero(ℝSAutF_N_star) : one(ℝSAutF_N_star) for i in 1:length(S), j in 1:length(S)]

    sos_pr_sym = LowCohomologySOS.sos_problem(
        M,
        order_unit,
        w_dec_matrix,
        length(collect(Σ)),
        1.0
    )

    # @info sos_pr_sym

    # TODO: more tests - especially the missing ones for sos_problem in a symmetrized version
# end