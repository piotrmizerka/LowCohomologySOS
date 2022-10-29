using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
using JuMP
using SymbolicWedderburn
using PermutationGroups
include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

function determine_transvection(g::Automorphism)
    @assert parent(g) isa Automorphism{<:FreeGroup}
    @assert length(word(g)) = 1

    A = alphabet(parent(g))

    return A[first(word(g))]
    end

Δ₁, Iₙ, half_basis, w_dec_matrix, Σ, basis, psd_basis, S = let half_radius = 2, N = 2
    SAut_F(n) = Groups.SpecialAutomorphismGroup(FreeGroup(n))
    SAut_F_N = SAut_F(N)

    S = let s = gens(SAut_F_N)
        [s; inv.(s)]
    end
    basis, sizes = Groups.wlmetric_ball(S, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]

    # @info length(basis)
    # @info length(half_basis)
    # @info length(S)

    F_SAut_F_2N = FreeGroup(length(S))

    quotient_hom = let source = F_SAut_F_2N, target = SAut_F_N
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([free_group_saut_index(i,length(S))]), source, target)
    end

    transvection_gen_dict = Dict([(determine_transvection(S[i]), gens(F_SAut_F_2N, i)) for i in 1:length(S)])
    ϱ(i,j, ε) = transvection_gen_dict[Groups.Transvection(:ϱ, i, j, ε)]
    λ(i,j, ε) = transvection_gen_dict[Groups.Transvection(:λ, i, j, ε)]
    ϱ(i,j) = ϱ(i,j,false)
    λ(i,j) = λ(i,j,false)

    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N
                     for j ∈ deleteat!(range_as_list, findall(j->j==i,range_as_list))]
    quadruples_1 = [(i,j,k,l) for k ∈ 1:N for l ∈ 1:N 
                              for i ∈ deleteat!(range_as_list, findall(i->i∈[k,l],range_as_list)) 
                              for j ∈ deleteat!(range_as_list, findall(j->j==k,range_as_list))]
    quadruples_2 = [(i,j,k,l) for k ∈ 1:N for l ∈ 1:N 
                              for i ∈ deleteat!(range_as_list, findall(i->i==l,range_as_list)) 
                              for j ∈ deleteat!(range_as_list, findall(j->j==k,range_as_list))]
    triples = [(i,j,k) for i ∈ 1:N
                       for j ∈ deleteat!(range_as_list, findall(j->j==i,range_as_list)) 
                       for k ∈ deleteat!(range_as_list, findall(k->k∈[i,j],range_as_list))]
    # TODO: I think these relation are in "first right, then left" composition convention - have to be changed!
    relations = vcat(
        [ϱ(i,j,ε)*ϱ(i,j,!ε) for (i,j) ∈ pairs for ε ∈ [true,false]],
        [λ(i,j,ε)*λ(i,j,ε) for (i,j) ∈ pairs for ε ∈ [true,false]],

        [ϱ(i,j,ε₁)*ϱ(k,l,ε₂)*ϱ(i,j,!ε₁)*ϱ(k,l,!ε₂) for (i,j,k,l) ∈ quadruples_1 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
        [λ(i,j,ε₁)*λ(k,l,ε₂)*λ(i,j,!ε₁)*λ(k,l,!ε₂) for (i,j,k,l) ∈ quadruples_1 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
        [ϱ(i,j,ε₁)*λ(k,l,ε₂)*ϱ(i,j,!ε₁)*λ(k,l,!ε₂) for (i,j,k,l) ∈ quadruples_2 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
        [λ(i,j,ε₁)*ϱ(k,l,ε₂)*λ(i,j,!ε₁)*ϱ(k,l,!ε₂) for (i,j,k,l) ∈ quadruples_2 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],

        [ϱ(i,j)*ϱ(j,k)*ϱ(i,j)^(-1)*ϱ(j,k)^(-1)*ϱ(i,k)^(-1) for (i,j,k) ∈ triples],
        [ϱ(i,j)*ϱ(j,k,true)*ϱ(i,j)^(-1)*ϱ(j,k,true)^(-1)*ϱ(i,k,true)^(-1) for (i,j,k) ∈ triples],
        [λ(i,j)*λ(j,k)*λ(i,j)^(-1)*λ(j,k)^(-1)*λ(i,k)^(-1) for (i,j,k) ∈ triples],
        [λ(i,j)*λ(j,k,true)*λ(i,j)^(-1)*λ(j,k,true)^(-1)*λ(i,k,true)^(-1) for (i,j,k) ∈ triples],
        [ϱ(i,j,true)*λ(j,k)*ϱ(i,j,true)^(-1)*λ(j,k)^(-1)*ϱ(i,k,true)^(-1) for (i,j,k) ∈ triples],
        [ϱ(i,j,true)*λ(j,k,true)*ϱ(i,j,true)^(-1)*λ(j,k,true)^(-1)*ϱ(i,k)^(-1) for (i,j,k) ∈ triples],
        [λ(i,j,true)*ϱ(j,k)*λ(i,j,true)^(-1)*ϱ(j,k)^(-1)*λ(i,k,true)^(-1) for (i,j,k) ∈ triples],
        [λ(i,j,true)*ϱ(j,k,true)*λ(i,j,true)^(-1)*ϱ(j,k,true)^(-1)*λ(i,k)^(-1) for (i,j,k) ∈ triples]
    )
    
    # @info length(relations)
    # @info "1"

    Δ₁, Iₙ = LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)

    # @info "2"

    Z_2_wr_S(n) = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(n))
    Σ = Z_2_wr_S(N)
    # w_dec_matrix = LowCohomologySOS.wedderburn_decomposition_matrix(Σ, constraints_basis, psd_basis, S)
    action = LowCohomologySOS.AlphabetPermutation(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)

    # @info "3"

    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)

    # @info length(w_dec_matrix.invariants)
    # @info w_dec_matrix.invariants[1]
    # @info w_dec_matrix.invariants[2345]

    Δ₁, Iₙ, half_basis, w_dec_matrix, Σ, basis, psd_basis, S
end

Δ₁_sos_problem = LowCohomologySOS.sos_problem_symmetrized(
    Δ₁, 
    Iₙ,
    w_dec_matrix,
    Σ,
    basis,
    psd_basis,
    S
)

SAut_F_4_data = (
    M = Δ₁,
    order_unit = Iₙ,
    half_basis = half_basis,
    RG = parent(first(Δ₁)),
)

solve_in_loop(
    Δ₁_sos_problem,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = SAut_F_4_data,
    w_dec_matrix
)
