using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
using SymbolicWedderburn
using PermutationGroups

include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

function group_data(half_radius, N, wreath_action)
    SAut_F(n) = Groups.SpecialAutomorphismGroup(FreeGroup(n))
    SAut_F_N = SAut_F(N)

    S_inv = let s = gens(SAut_F_N)
        [s; inv.(s)]
    end
    S = (wreath_action ? S_inv : gens(SAut_F_N))
    basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]
    ℝSAutF_N_star = LowCohomologySOS.group_ring(SAut_F_N, half_basis, star_multiplication = true)

    return SAut_F_N, ℝSAutF_N_star.basis, half_basis, S
end

function wedderburn_data(basis, half_basis, S, N, wreath_action)
    @time begin
        if wreath_action
            Z_2_wr_S(n) = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(n))
            Σ = Z_2_wr_S(N)
        else
            Σ = PermutationGroups.SymmetricGroup(N)
        end
        actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)
        constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    end

    return constraints_basis, psd_basis, Σ, actions
end

function determine_transvection(g)
    @assert length(word(g)) == 1
    
    A = alphabet(parent(g))

    return A[first(word(g))]
end

function free_group_saut_index(i::Integer, F_G, S)
    gen_id = (-1, false)
    for j in eachindex(gens(F_G))
        if gens(F_G, j) == F_G([i])
            gen_id = (j, false)
        elseif gens(F_G, j) == inv(F_G([i]))
            gen_id = (j, true)
        end
    end
    
    return gen_id[2] ? word(inv(S[gen_id[1]]))[1] : word(S[gen_id[1]])[1]
end

const half_radius = 2;
const N = 4;
const wreath_action = false;

SAut_F_N, basis, half_basis, S = group_data(half_radius, N, wreath_action)

Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = let
    if !wreath_action
        F_G = FreeGroup(alphabet(SAut_F_N))
        quotient_hom = let source = F_G, target = SAut_F_N
            Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
        end
    else
        F_G = FreeGroup(length(S))
        quotient_hom = let source = F_G, target = SAut_F_N
            Groups.Homomorphism((i, F, G) -> Groups.word_type(G)(
                [free_group_saut_index(i,F_G, S)]), source, target)
        end
    end

    # check if the quotient homomorphism is defined properly
    @assert length(gens(F_G)) == length(S)
    for i in 1:length(S)
        @assert quotient_hom(gens(F_G,i)) == S[i]
    end

    transvection_gen_dict = Dict([(determine_transvection(S[i]), gens(F_G, i)) for i in 1:length(S)])
    ϱ(i,j, ε) = transvection_gen_dict[Groups.Transvection(:ϱ, i, j, ε)]
    λ(i,j, ε) = transvection_gen_dict[Groups.Transvection(:λ, i, j, ε)]

    range_as_list = [i for i in 1:N]
    quadruples_1 = [(i,j,k,l) for k ∈ 1:N
                            for l ∈ deleteat!(copy(range_as_list), findall(l->l==k,copy(range_as_list))) 
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i∈[k,l],copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    quadruples_2 = [(i,j,k,l) for k ∈ 1:N
                            for l ∈ deleteat!(copy(range_as_list), findall(l->l==k,copy(range_as_list)))
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i==l,copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ 1:N
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]

    # The relations are derived from the Gersten's article, https://www.sciencedirect.com/science/article/pii/0022404984900628,
    # Theorem 2.8. Note that our convention assumes the automorphism composition order reversed with respect to Gersten's.
    # Therefore, the order of letters in the relators had to be reversed as well (earlier, we had to change Gertsen's symbols E_a_b
    # to the ϱ and λ notation used in https://annals.math.princeton.edu/2021/193-2/p03).
    if wreath_action
        pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]

        relations = vcat(
            [ϱ(i,j,ε)*ϱ(i,j,!ε) for (i,j) ∈ pairs for ε ∈ [true,false]],
            [λ(i,j,ε)*λ(i,j,!ε) for (i,j) ∈ pairs for ε ∈ [true,false]],

            [ϱ(k,l,!ε₂)*ϱ(i,j,!ε₁)*ϱ(k,l,ε₂)*ϱ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_1 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
            [λ(k,l,!ε₂)*λ(i,j,!ε₁)*λ(k,l,ε₂)*λ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_1 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
            [ϱ(k,l,!ε₂)*λ(i,j,!ε₁)*ϱ(k,l,ε₂)*λ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_2 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
            [λ(k,l,!ε₂)*ϱ(i,j,!ε₁)*λ(k,l,ε₂)*ϱ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_2 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],

            [ϱ(i,k)*ϱ(j,k)*ϱ(i,j,true)*ϱ(j,k,true)*ϱ(i,j) for (i,j,k) ∈ triples],
            [ϱ(i,k,true)*ϱ(j,k,true)*ϱ(i,j,true)*ϱ(j,k)*ϱ(i,j) for (i,j,k) ∈ triples],
            [λ(i,k)*λ(j,k)*λ(i,j,true)*λ(j,k,true)*λ(i,j) for (i,j,k) ∈ triples],
            [λ(i,k,true)*λ(j,k,true)*λ(i,j,true)*λ(j,k)*λ(i,j) for (i,j,k) ∈ triples],
            [ϱ(i,k)*λ(j,k,true)*ϱ(i,j)*λ(j,k)*ϱ(i,j,true) for (i,j,k) ∈ triples],
            [ϱ(i,k,true)*λ(j,k)*ϱ(i,j)*λ(j,k,true)*ϱ(i,j,true) for (i,j,k) ∈ triples],
            [λ(i,k)*ϱ(j,k,true)*λ(i,j)*ϱ(j,k)*λ(i,j,true) for (i,j,k) ∈ triples],
            [λ(i,k,true)*ϱ(j,k)*λ(i,j)*ϱ(j,k,true)*λ(i,j,true) for (i,j,k) ∈ triples]
        )
    else
        ϱ(i,j) = ϱ(i,j,false)
        λ(i,j) = λ(i,j,false)

        relations = vcat(
            [ϱ(k,l)^(-1)*ϱ(i,j)^(-1)*ϱ(k,l)*ϱ(i,j) for (i,j,k,l) ∈ quadruples_1],
            [λ(k,l)^(-1)*λ(i,j)^(-1)*λ(k,l)*λ(i,j) for (i,j,k,l) ∈ quadruples_1],
            [λ(k,l)^(-1)*ϱ(i,j)^(-1)*λ(k,l)*ϱ(i,j) for (i,j,k,l) ∈ quadruples_2],

            [ϱ(i,k)^(-1)*ϱ(j,k)^(-1)*ϱ(i,j)^(-1)*ϱ(j,k)*ϱ(i,j) for (i,j,k) ∈ triples],
            [ϱ(i,k)*ϱ(j,k)*ϱ(i,j)^(-1)*ϱ(j,k)^(-1)*ϱ(i,j) for (i,j,k) ∈ triples],
            [λ(i,k)^(-1)*λ(j,k)^(-1)*λ(i,j)^(-1)*λ(j,k)*λ(i,j) for (i,j,k) ∈ triples],
            [λ(i,k)*λ(j,k)*λ(i,j)^(-1)*λ(j,k)^(-1)*λ(i,j) for (i,j,k) ∈ triples],
            [ϱ(i,k)*λ(j,k)^(-1)*ϱ(i,j)*λ(j,k)*ϱ(i,j)^(-1) for (i,j,k) ∈ triples],
            [ϱ(i,k)^(-1)*λ(j,k)*ϱ(i,j)*λ(j,k)^(-1)*ϱ(i,j)^(-1) for (i,j,k) ∈ triples],
            [λ(i,k)*ϱ(j,k)^(-1)*λ(i,j)*ϱ(j,k)*λ(i,j)^(-1) for (i,j,k) ∈ triples],
            [λ(i,k)^(-1)*ϱ(j,k)*λ(i,j)*ϱ(j,k)^(-1)*λ(i,j)^(-1) for (i,j,k) ∈ triples]
        )
    end

    LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)
end

constraints_basis, psd_basis, Σ, action = wedderburn_data(basis, half_basis, S, N, wreath_action);

@time "\tWedderburn total" begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    Δ₁_sos_problem = LowCohomologySOS.sos_problem(
        Δ₁, 
        Iₙ,
        w_dec_matrix,
        length(collect(Σ)),
        1.0
    )
end

SAut_F_N_data = (
    M = Δ₁,
    order_unit = Iₙ,
    half_basis = half_basis,
    RG = parent(first(Δ₁)),
)

solve_in_loop(
    Δ₁_sos_problem,
    w_dec_matrix,
    logdir = "./LowCohomologySOS/logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 10_000),
    data = SAut_F_N_data
)
