# We want to find alpha and beta such that
# αd₁^*d₁+βd₀d₀^* = symmetrized induced Laplacians from SL(3,Z) to SL(4,Z)
using Revise 

using StarAlgebras
using Groups
using LowCohomologySOS

const N = 3
const M = 4

i = LowCohomologySOS.sln_slm_embedding(N, M)

slN = i.source
slM = i.target

const half_radius = 2

slN_S = gens(slN)
slN_S_inv = let s = slN_S
    [s; inv.(s)]
end
slM_S = gens(slM)
slM_S_inv = let s = slM_S
    [s; inv.(s)]
end

slN_half_basis, slN_sizes = Groups.wlmetric_ball(slN_S_inv, radius = half_radius)
slM_half_basis, slM_sizes = Groups.wlmetric_ball(slM_S_inv, radius = half_radius)

slN_Δ₁, slN_Iₙ, slN_Δ₁⁺, slN_Δ₁⁻ = LowCohomologySOS.sln_laplacians(slN, slN_half_basis, slN_S)
slM_Δ₁, slM_Iₙ, slM_Δ₁⁺, slM_Δ₁⁻ = LowCohomologySOS.sln_laplacians(slM, slM_half_basis, slM_S)

@assert parent(first(slM_Δ₁⁺)) == parent(first(slM_Δ₁⁻))

RG_prime = parent(first(slM_Δ₁⁺))
Δ₁_emb = LowCohomologySOS.embed_matrix(slN_Δ₁, i, RG_prime)

@assert parent(first(Δ₁_emb)) == parent(first(slM_Δ₁⁻))

using PermutationGroups

Δ₁_emb_symmetrized = let # n = 4
    # Σ = PermutationGroups.SymmetricGroup(n)
    Σ = PermGroup(Perm{Int8}[perm"(1,2,3)", perm"(1,2)(3,4)"]) # alternating group A₄
    LowCohomologySOS.weyl_symmetrize_matrix(Δ₁_emb, Σ, LowCohomologySOS._conj)
end

slM_Δ₁⁺[1,1]
slM_Δ₁⁻[1,1]
Δ₁_emb_symmetrized[1,1]

M_ = 3*slM_Δ₁⁺+3*slM_Δ₁⁻-Δ₁_emb_symmetrized

S = gens(slM)
S_inv = let s = S
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 1)
RG_prime_reduced = LowCohomologySOS.group_ring(slM, half_basis, star_multiplication = true)
function id__(letter_id, SLₙℤ, G)
    return word(SLₙℤ([letter_id]))
end
id_ = let source = slM, target = slM
    Groups.Homomorphism(id__, source, target, check = false)
end

M__ = LowCohomologySOS.embed_matrix(M_, id_, RG_prime_reduced)
slM_Iₙ__ = LowCohomologySOS.embed_matrix(slM_Iₙ, id_, RG_prime_reduced)

_data = (
    M = M__,
    order_unit = slM_Iₙ__,
    half_basis = half_basis,
    RG = parent(first(M__)),
)

M_sos_problem = LowCohomologySOS.sos_problem(M__, slM_Iₙ__)

solve_in_loop(
    M_sos_problem,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = _data
)

# the code below may be not that helpful in fact ... ###################################################
using JuMP

function alpha_beta_problem(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    M::AbstractMatrix{<:AlgebraElement},
    alpha_ge_beta = true
)
    @assert size(A) == size(B) == size(M)
    @assert !isempty(A) && !isempty(B) && !isempty(M)

    RG = parent(first(M))

    @assert all(x -> parent(x) === RG, A)
    @assert all(x -> parent(x) === RG, B)
    @assert all(x -> parent(x) === RG, M)

    result = JuMP.Model()

    JuMP.@variable(result, α)
    JuMP.@variable(result, β)

    # The perfect situation would be in the case α = β
    if alpha_ge_beta 
        JuMP.@constraint(result, α >= β)
        JuMP.@objective(result, Max, β-α)
    else
        JuMP.@constraint(result, β >= α)
        JuMP.@objective(result, Max, α-β)
    end

    for idx in CartesianIndices(M)
        aij = A[idx]
        bij = B[idx]
        mij = M[idx]

        for g in basis(RG)
            JuMP.@constraint(result, α *aij(g)+β*bij(g) == mij(g))
        end
    end

    return result
end

include(joinpath(@__DIR__, "optimizers.jl"))

function alpha_beta_adjust(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    M::AbstractMatrix{<:AlgebraElement},
    alpha_ge_beta = true,
    optimizer = scs_opt(eps = 1e-9, max_iters = 10_000)
)
    problem = alpha_beta_problem(A, B, M, alpha_ge_beta)

    JuMP.set_optimizer(problem, optimizer)
    JuMP.optimize!(problem)

    return JuMP.value.(problem[:α]), JuMP.value.(problem[:β])
end

# Code below intended for tests ###############################################################
function cyclic_group(n::Integer)
    A = Alphabet([:a, :A], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e])

    return Cₙ
end

G = cyclic_group(3)

using LowCohomologySOS

RG =  LowCohomologySOS.group_ring(G, 1)

a, = Groups.gens(G)
A = [one(RG) RG(a); one(RG) one(RG)]
B = [one(RG) zero(RG); zero(RG) one(RG)]
Mx= [2*one(RG) RG(a); one(RG) 2*one(RG)]

alpha_beta_adjustx = alpha_beta_adjust(A, B, Mx)
alpha_beta_adjustx = alpha_beta_adjust(A, B, Mx, false)
#############################################################################################
# The last functions to call:
alpha_beta_adjust(slM_Δ₁⁺, slM_Δ₁⁻, Δ₁_emb_symmetrized)
alpha_beta_adjust(slM_Δ₁⁺, slM_Δ₁⁻, Δ₁_emb_symmetrized, false)
###################################################################################################