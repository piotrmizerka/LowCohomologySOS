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
            eσ_inv, fσ_inv = LowCohomologySOS.action_on_group(action, σ^(-1), e), LowCohomologySOS.action_on_group(action, σ^(-1), f)
            id_eσ_inv, id_fσ_inv = half_basis_idies[eσ_inv], half_basis_idies[fσ_inv]
            result[id_e, id_f] = P[id_eσ_inv, id_fσ_inv]
        end
    end

    return result
end

const N = 2
const half_radius = 1

S, ℝSAutF_N_star, half_basis, w_dec_matrix, A_gs_cart, Σ, psd_basis_length, actions = let
    SAutF_N = Groups.SpecialAutomorphismGroup(FreeGroup(N))

    S = let s = Groups.gens(SAutF_N)
        [s; inv.(s)]
    end
    S = unique!(S)

    basis, sizes = Groups.wlmetric_ball(S, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]
    ℝSAutF_N_star = LowCohomologySOS.group_ring(SAutF_N, half_basis, star_multiplication = true)

    Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(N))
    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, ℝSAutF_N_star.basis)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(ℝSAutF_N_star.basis, half_basis, S)
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, actions, constraints_basis, psd_basis)

    cnstrs =  dropzeros!.(sparse.(LowCohomologySOS.constraints(ℝSAutF_N_star.mstructure)))
    A_gs_cart = [findall(!iszero,c) for c in cnstrs]

    S, ℝSAutF_N_star, half_basis, w_dec_matrix, A_gs_cart, Σ, length(psd_basis), actions
end

@testset "invariant_constraint_matrix" begin
    inv_cnstr_matrix = LowCohomologySOS.invariant_constraint_matrix(
        rand(w_dec_matrix.invariants),
        A_gs_cart,
        length(half_basis),
        length(collect(Σ))
    )

    @test size(inv_cnstr_matrix) == (psd_basis_length, psd_basis_length)

    half_basis_idies = Dict(half_basis[i] => i for i in 1:length(half_basis))
    S_idies = Dict(S[i] => i for i in 1:length(S))
    for i in 1:length(S)
        for j in 1:length(S)
            lhs = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(i, j)]
            for σ ∈ Σ
                pbe_σ_inv = SymbolicWedderburn.action(actions.alphabet_perm, σ^(-1), LowCohomologySOS.PSDBasisElement(S[i], S[j]))
                iσ_inv, jσ_inv = S_idies[pbe_σ_inv.s], S_idies[pbe_σ_inv.g]
                P = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(iσ_inv, jσ_inv)]
                rhs = act_on_matrix(P, σ, half_basis, half_basis_idies, actions.alphabet_perm)
                @test lhs == rhs
            end
        end
    end
end
 
@testset "symmetrized matrix SOS problem" begin
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
end