using SparseArrays
using Kronecker
using LinearAlgebra

# act on real valued matrix defined by a matrix over a group ring
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

function wedderburn_tests(action_type, group_name, N, half_radius)
    if group_name == "sautfn"
        G = Groups.SpecialAutomorphismGroup(FreeGroup(N))
    elseif group_name == "sln"
        G = MatrixGroups.SpecialLinearGroup{N}(Int8)
    end

    S_inv = let S = gens(G)
        [S; inv.(S)]
    end

    if action_type == "wreath"
        S = S_inv
        Σ = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(N))
    elseif action_type == "symmetric"
        S = gens(G)
        Σ = PermutationGroups.SymmetricGroup(N)
    end

    basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]
    RG_star = LowCohomologySOS.group_ring(G, half_basis, star_multiplication = true)

    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, RG_star.basis)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(RG_star.basis, half_basis, S)
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, actions, constraints_basis, psd_basis)

    cnstrs =  dropzeros!.(sparse.(LowCohomologySOS.constraints(RG_star.mstructure)))
    A_gs_cart = [findall(!iszero,c) for c in cnstrs]

    invariant_testset_name = action_type*" invariant_constraint_matrix for "*group_name
    @testset "$invariant_testset_name" begin
        inv_cnstr_matrix = LowCohomologySOS.invariant_constraint_matrix(
            rand(w_dec_matrix.invariants),
            A_gs_cart,
            length(half_basis)
        )

        @test size(inv_cnstr_matrix) == (length(psd_basis), length(psd_basis))

        half_basis_idies = Dict(half_basis[i] => i for i in 1:length(half_basis))
        S_idies = Dict(S[i] => i for i in 1:length(S))
        for i in 1:length(S)
            for j in 1:length(S)
                lhs = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(i, j)]
                for σ ∈ Σ
                    pbe_σ_inv = SymbolicWedderburn.action(actions.alphabet_perm, σ^(-1), LowCohomologySOS.PSDBasisElement(S[i], S[j]))
                    iσ_inv, jσ_inv = S_idies[pbe_σ_inv.generator], S_idies[pbe_σ_inv.basis_elt]
                    P = @view inv_cnstr_matrix[LowCohomologySOS.KroneckerDelta{length(S)}(iσ_inv, jσ_inv)]
                    rhs = act_on_matrix(P, σ, half_basis, half_basis_idies, actions.alphabet_perm)
                    @test lhs == rhs
                end
            end
        end
    end

    sos_problem_testset_name = action_type*" matrix SOS problem for "*group_name
    @testset "$sos_problem_testset_name" begin
        function test_solution(M)
            # there is no point of finding a solution if we don't provide invariant matrix
            for σ in Σ
                @assert LowCohomologySOS.act_on_matrix(M, σ, actions.alphabet_perm, S) == M
            end

            sos_pr_sym = LowCohomologySOS.sos_problem(
                M,
                order_unit,
                w_dec_matrix,
                1.0
            )
        
            JuMP.set_optimizer(sos_pr_sym[1], scs_opt(eps = 1e-7, max_iters = 5_000))
            JuMP.optimize!(sos_pr_sym[1])
        
            λ, Q = LowCohomologySOS.get_solution(sos_pr_sym[1], sos_pr_sym[2], w_dec_matrix)
        
            estimated_sos = LowCohomologySOS.sos_from_matrix(RG_star, Q, half_basis)
            eoi = M-λ*order_unit
            residual = eoi - estimated_sos
            l1_norm = sum(x -> norm(x, 1), residual)
        
            @test l1_norm.hi < 0.001
        end

        order_unit = [i ≠ j ? zero(RG_star) : one(RG_star) for i in 1:length(S), j in 1:length(S)]

        M = 2*order_unit # check if it finds a solution for identity - it should at all times :)
        test_solution(M)

        M = [one(RG_star) for i in 1:length(S), j in 1:length(S)]+100*order_unit
        test_solution(M)

        invariant_elt = sum(RG_star(s) for s in S)

        ξ = [i ≠ j ? zero(RG_star) : one(RG_star)+invariant_elt for i in 1:length(S), j in 1:length(S)]
        M = ξ^2+0.7*order_unit
        test_solution(M)

        ξ = [i ≠ j ? zero(RG_star) : one(RG_star)+invariant_elt for i in 1:length(S), j in 1:length(S)]
        M = ξ^2+0.7*order_unit
        test_solution(M)
    end
end

function sq_test(action_type)
    testset_name = action_type*" Sq cetification for SL₃(ℤ)"
    @testset "$testset_name" begin 
        N = 3
        G = MatrixGroups.SpecialLinearGroup{N}(Int8)
        S_inv = let S = gens(G)
            [S; inv.(S)]
        end
        S = (action_type == "symmetric") ? gens(G) : S_inv
        basis, sizes = Groups.wlmetric_ball(S_inv, radius = 4)
        half_basis = basis[1:sizes[2]]
        Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(G, half_basis, S)
        half_basis_wed = basis[1:sizes[1]]
        RG_wed = LowCohomologySOS.group_ring(G, half_basis_wed, star_multiplication = true)

        Δ₁⁻_wed = LowCohomologySOS.embed.(identity, Δ₁⁻, Ref(RG_wed))
        sq, adj, op = LowCohomologySOS.sq_adj_op(Δ₁⁻_wed, S)
        Iₙ = [i ≠ j ? zero(RG_wed) : one(RG_wed) for i in 1:length(S), j in 1:length(S)]

        Σ = (action_type == "symmetric") ? PermutationGroups.SymmetricGroup(N) : Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(N))
        actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, RG_wed.basis)
        constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(RG_wed.basis, half_basis_wed, S)
        w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, actions, constraints_basis, psd_basis)

        sos_pr_sym = LowCohomologySOS.sos_problem(
            sq,
            Iₙ,
            w_dec_matrix,
            1.0
        )
        JuMP.set_optimizer(sos_pr_sym[1], scs_opt(eps = 1e-7, max_iters = 5_000))
        JuMP.optimize!(sos_pr_sym[1])
        λ, Q = LowCohomologySOS.get_solution(sos_pr_sym[1], sos_pr_sym[2], w_dec_matrix)
        estimated_sos = LowCohomologySOS.sos_from_matrix(RG_wed, Q, half_basis_wed)
        eoi = sq-λ*Iₙ
        residual = eoi - estimated_sos
        l1_norm = sum(x -> norm(x, 1), residual)
    
        @test l1_norm.hi < 0.001
    end
end

wedderburn_tests("symmetric", "sln", 3, 1)
wedderburn_tests("symmetric", "sautfn", 3, 1)
wedderburn_tests("wreath", "sln", 3, 1)
wedderburn_tests("wreath", "sautfn", 3, 1)
sq_test("symmetric")
sq_test("wreath")
