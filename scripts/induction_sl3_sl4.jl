# We want to find alpha and beta such that
# αd₁^*d₁+βd₀d₀^* = symmetrized induced Laplacians from SL(3,Z) to SL(4,Z)
using Revise 

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using StarAlgebras
using Groups
using LowCohomologySOS

const N = 3
const M = 4

SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)

slN, slM = SL(N,Int8), SL(M,Int8)

const half_radius = 2

slN_S = gens(slN);
slN_S_inv = let s = slN_S
    [s; inv.(s)]
end;
slM_S = gens(slM);
slM_S_inv = let s = slM_S
    [s; inv.(s)]
end

slN_half_basis, slN_sizes = Groups.wlmetric_ball(slN_S_inv, radius = half_radius)
slM_half_basis, slM_sizes = Groups.wlmetric_ball(slM_S_inv, radius = half_radius)

slN_Δ₁, slN_Iₙ, slN_Δ₁⁺, slN_Δ₁⁻ = LowCohomologySOS.sln_laplacians(slN, slN_half_basis, slN_S);
slM_Δ₁, slM_Iₙ, slM_Δ₁⁺, slM_Δ₁⁻ = LowCohomologySOS.sln_laplacians(slM, slM_half_basis, slM_S);

using PermutationGroups

# Σ = PermutationGroups.SymmetricGroup(M)
# Σ = PermGroup(Perm{Int8}[perm"(1,3,5)", perm"(1,2)(3,4)"]) # alternating group A₅
Σ = PermGroup(Perm{Int8}[perm"(1,2,3)", perm"(1,2)(3,4)"]) # alternating group A₄

RG = parent(first(slN_Δ₁⁺))
RG_prime = parent(first(slM_Δ₁⁺))

I = [i ≠ j ? zero(RG) : one(RG) for i in eachindex(gens(slN)), j in eachindex(gens(slN))]
I_emb_symmetrized = LowCohomologySOS.embedded_symmetrized_matrix(I, Σ, RG_prime)
ones = [one(RG) for i in eachindex(gens(slN)), j in eachindex(gens(slN))]
ones_emb_symmetrized = LowCohomologySOS.embedded_symmetrized_matrix(ones, Σ, RG_prime)

elts = collect(Σ)
elts[1]
elts[2]
elts[3]
elts[4]
elts[5]
elts[6]
elts[7]
elts[8]
elts[9]
elts[10]
elts[11]
elts[12]

using Combinatorics

subsets = collect(powerset(elts))
length(subsets)

good_subset = [-1]
it = 1
for subset in subsets
    if it % 1000 == 0
        println(it)
    end
    it += 1
    if length(subset) > 3
        result = [zero(RG_prime) for i in eachindex(gens(slM)), j in eachindex(gens(slM))]
        for σ in subset
            iσ = LowCohomologySOS.sln_slm_embedding(slN, slM, σ)
            result += LowCohomologySOS.embed_matrix(ones, iσ, RG_prime);
        end
        coeffs = []
        for i in 1:12
            for j in 1:12
                append!(coeffs,result[i,j](one(slM)))
            end
        end
        if length(Set(coeffs)) == 2
            good_subset = subset
            break
        end
    end
end

good_subset




show(IOContext(stdout, :limit => true, :displaysize => (120, 120)), "text/plain", result)

Δ₁⁺_emb_symmetrized = LowCohomologySOS.embedded_symmetrized_matrix(slN_Δ₁⁺, Σ, RG_prime)
Δ₁⁻_emb_symmetrized = LowCohomologySOS.embedded_symmetrized_matrix(slN_Δ₁⁻, Σ, RG_prime)

@assert parent(first(Δ₁⁺_emb_symmetrized)) == parent(first(Δ₁⁻_emb_symmetrized)) == parent(first(slM_Δ₁⁺)) == parent(first(slM_Δ₁⁻))

slM_Δ₁⁻
Δ₁⁻_emb_symmetrized
72*slM_Δ₁⁻-Δ₁⁻_emb_symmetrized

using JuMP

function alpha_beta_gamma_problem(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    C::AbstractMatrix{<:AlgebraElement},
    D::AbstractMatrix{<:AlgebraElement}
)
    @assert size(A) == size(B) == size(C) == size(D)
    @assert !isempty(A) && !isempty(B) && !isempty(C) && !isempty(D)

    RG = parent(first(A))

    @assert all(x -> parent(x) === RG, A)
    @assert all(x -> parent(x) === RG, B)
    @assert all(x -> parent(x) === RG, C)
    @assert all(x -> parent(x) === RG, D)

    result = JuMP.Model()

    JuMP.@variable(result, α)
    JuMP.@variable(result, β)
    JuMP.@variable(result, γ)

    # JuMP.@constraint(result, α >= 0.01)
    # JuMP.@constraint(result, β >= 0.01)
    # JuMP.@constraint(result, γ >= 0.01)

    for idx in CartesianIndices(A)
        aij = A[idx]
        bij = B[idx]
        cij = C[idx]
        dij = D[idx]

        for g in basis(RG)
            JuMP.@constraint(result, α*aij(g) + β*bij(g) == γ*cij(g) + dij(g))
        end
    end

    return result
end

include(joinpath(@__DIR__, "optimizers.jl"))

function alpha_beta_adjust(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    C::AbstractMatrix{<:AlgebraElement},
    D::AbstractMatrix{<:AlgebraElement},
    optimizer = scs_opt(eps = 1e-9, max_iters = 10_000)
)
    problem = alpha_beta_gamma_problem(A, B, C, D)

    JuMP.set_optimizer(problem, optimizer)
    JuMP.optimize!(problem)

    return JuMP.value.(problem[:α]), JuMP.value.(problem[:β]), JuMP.value.(problem[:γ])
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
C = [2*one(RG) RG(a); one(RG) 2*one(RG)]
D = [2*one(RG) RG(a); one(RG) 2*one(RG)]

alpha_beta_adjustx = alpha_beta_adjust(A, B, C, D)
#############################################################################################

# The last functions to call:
alpha_beta_adjust(slM_Δ₁⁺, slM_Δ₁⁻, Δ₁⁺_emb_symmetrized, Δ₁⁻_emb_symmetrized)
