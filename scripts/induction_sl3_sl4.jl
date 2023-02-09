# We want to find alpha and beta such that
# αd₁^*d₁+βd₀d₀^* = symmetrized induced Laplacians from SL(3,Z) to SL(4,Z)
using Revise 

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using StarAlgebras
using Groups
using LowCohomologySOS

const N = 4
const M = 5

i = LowCohomologySOS.sln_slm_embedding(N, M)

slN = i.source
slM = i.target

const half_radius = 2

slN_S = gens(slN);
slN_S_inv = let s = slN_S
    [s; inv.(s)]
end;
slM_S = gens(slM);
slM_S_inv = let s = slM_S
    [s; inv.(s)]
end;

slN_half_basis, slN_sizes = Groups.wlmetric_ball(slN_S_inv, radius = half_radius);
slM_half_basis, slM_sizes = Groups.wlmetric_ball(slM_S_inv, radius = half_radius);

slN_Δ₁, slN_Iₙ, slN_Δ₁⁺, slN_Δ₁⁻ = LowCohomologySOS.laplacians(slN, slN_half_basis, slN_S, sq_adj_op_ = "adj");
slM_Δ₁, slM_Iₙ, slM_Δ₁⁺, slM_Δ₁⁻ = LowCohomologySOS.laplacians(slM, slM_half_basis, slM_S, sq_adj_op_ = "adj");
slN_sq, slN_adj, slN_op = LowCohomologySOS.sq_adj_op(slN_Δ₁⁻, slN_S)
slM_sq, slM_adj, slM_op = LowCohomologySOS.sq_adj_op(slM_Δ₁⁻, slM_S)

RG_prime = parent(first(slM_Δ₁⁺))

Δ₁⁺_emb = LowCohomologySOS.embed_matrix(slN_Δ₁⁺, i, RG_prime)
adj_emb = LowCohomologySOS.embed_matrix(slN_adj, i, RG_prime)

@assert parent(first(Δ₁⁺_emb)) == parent(first(adj_emb)) == parent(first(slM_Δ₁⁺)) == parent(first(slM_adj))

using PermutationGroups

Δ₁⁺_emb_symmetrized = let
    # Σ = PermutationGroups.SymmetricGroup(4)
    # Σ = PermGroup(Perm{Int8}[perm"(1,2,3)", perm"(1,2)(3,4)"]) # alternating group A₄
    Σ = PermGroup(Perm{Int8}[perm"(1,2,3)", perm"(1,2,3,4,5)"]) # alternating group A₅
    LowCohomologySOS.weyl_symmetrize_matrix(Δ₁⁺_emb, Σ, LowCohomologySOS._conj, slM_S)
end

adj_emb_symmetrized = let
    # Σ = PermutationGroups.SymmetricGroup(4)
    # Σ = PermGroup(Perm{Int8}[perm"(1,2,3)", perm"(1,2)(3,4)"]) # alternating group A₄
    Σ = PermGroup(Perm{Int8}[perm"(1,2,3)", perm"(1,2,3,4,5)"]) # alternating group A₅
    LowCohomologySOS.weyl_symmetrize_matrix(adj_emb, Σ, LowCohomologySOS._conj, slM_S)
end

24*slM_Δ₁⁺-Δ₁⁺_emb_symmetrized # it looks like symmetrization works for upper Laplacians!
24*slM_adj-adj_emb_symmetrized # Adj symmetrizes as well with the same pace!!
