# We want to find alpha and beta such that
# αd₁^*d₁+βd₀d₀^* = symmetrized induced Laplacians from SL(3,Z) to SL(4,Z)
using Revise 

using StarAlgebras
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

using Groups

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
M = [2*one(RG) RG(a); one(RG) 2*one(RG)]

alpha_beta_adjustx = alpha_beta_adjust(A, B, M)
alpha_beta_adjustx = alpha_beta_adjust(A, B, M, false)
#############################################################################################

using LowCohomologySOS

const N = 4

i = LowCohomologySOS.sln_slm_embedding(3, N)

slN = i.target

const half_radius = 2

Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.sln_laplacians(slN, half_radius)

function laplacian_embedding(n::Integer, m::Integer)
    # TODO ??
end

# Laplacian embedding from SL(3,Z) to SL(4,Z) ###############################################
@assert parent(first(Δ₁⁺)) == parent(first(Δ₁⁻))

RG_prime = parent(first(Δ₁⁺))

Δ₁_emb = let half_radius = 2
    sl3 = i.source
    S = let s = gens(sl3)
        [s; inv.(s)]
    end
    half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)

    F_sl_3_z = FreeGroup(alphabet(sl3))
    e12, e13, e21, e23, e31, e32 = Groups.gens(F_sl_3_z)

    quotient_hom = let source = F_sl_3_z, target = sl3
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    relations = [
        e12 * e13 * e12^(-1) * e13^(-1),
        e12 * e32 * e12^(-1) * e32^(-1),
        e13 * e23 * e13^(-1) * e23^(-1),
        e23 * e21 * e23^(-1) * e21^(-1),
        e21 * e31 * e21^(-1) * e31^(-1),
        e31 * e32 * e31^(-1) * e32^(-1),
        e12 * e23 * e12^(-1) * e23^(-1) * e13^(-1),
        e13 * e32 * e13^(-1) * e32^(-1) * e12^(-1),
        e21 * e13 * e21^(-1) * e13^(-1) * e23^(-1),
        e23 * e31 * e23^(-1) * e31^(-1) * e21^(-1),
        e31 * e12 * e31^(-1) * e12^(-1) * e32^(-1),
        e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1)
    ]

    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)

    Δ₁_emb = LowCohomologySOS.embed_matrix(Δ₁, i, RG_prime)
    Δ₁_emb
end

@assert parent(first(Δ₁_emb)) == parent(first(Δ₁⁻))
#############################################################################################

using SymbolicWedderburn
using SparseArrays

function act_on_matrix(
    M::AbstractMatrix{<:AlgebraElement}, # M has to be square and indexed by the generators
    σ::Groups.GroupElement,
    act::LowCohomologySOS.AlphabetPermutation
)
    RG = parent(first(M))
    G = parent(first(basis(RG)))
    S = gens(G)
    gen_idies = Dict(S[i] => i for i in eachindex(S))
    basis_ = basis(RG)

    result = [zero(RG) for i in eachindex(S), j in eachindex(S)]

    for i in eachindex(S)
        for j in eachindex(S)
            s, t = S[i], S[j]
            s_σ_inv, t_σ_inv = SymbolicWedderburn.action(act, σ^(-1), s), SymbolicWedderburn.action(act, σ^(-1), t)
            # s_σ_inv, t_σ_inv = SymbolicWedderburn.action(act, σ, s), SymbolicWedderburn.action(act, σ, t) # for left action
            coeffs_ = coeffs(M[gen_idies[s_σ_inv],gen_idies[t_σ_inv]])
            for ind in SparseArrays.nonzeroinds(coeffs_)
                result[i,j] += coeffs_[ind]*RG(SymbolicWedderburn.action(act, σ, basis_[ind]))
                # result[i,j] += coeffs_[ind]*RG(SymbolicWedderburn.action(act, σ^(-1), basis_[ind])) # for left action
            end
        end
    end

    return result
end

function alphabet_permutation(
    A::Alphabet, 
    G, 
    op,
    n::Integer
)
    return LowCohomologySOS.AlphabetPermutation(
        Dict(
            g => PermutationGroups.Perm([A[op(l, g, n)] for l in A.letters]) for
            g in G
        ),
    )
end

using PermutationGroups

function sln_alphabet_op(
    l,
    σ::PermutationGroups.AbstractPerm,
    n::Integer
)
    return Groups.MatrixGroups.ElementaryMatrix{n}(l.i^σ, l.j^σ, l.val)
end

# Code below intended for tests (change @assert to @test) ###############################################################
const N = 2

sln = SL(N,Int8)
Σ = PermutationGroups.SymmetricGroup(N)
alphabet_permutation_ = alphabet_permutation(alphabet(sln), Σ, sln_alphabet_op, N)

RG =  LowCohomologySOS.group_ring(sln, 1)
e12, e21 = gens(sln)
M = [RG(e12) zero(RG); RG(e21) one(RG)]

σ = collect(Σ)[2]

m_σ = act_on_matrix(M, σ, alphabet_permutation_)

# action definition agreement:
const n = 3

sln = SL(n,Int8)
Σ = PermutationGroups.SymmetricGroup(n)
alphabet_permutation_ = alphabet_permutation(alphabet(sln), Σ, sln_alphabet_op, n)

RG =  LowCohomologySOS.group_ring(sln, 1)
e12, e13, e21, e23, e31, e32 = gens(sln)
M = [
    RG(e12) zero(RG) zero(RG) zero(RG) zero(RG) zero(RG);
    zero(RG) zero(RG) zero(RG) zero(RG) zero(RG) zero(RG); 
    zero(RG) zero(RG) zero(RG) zero(RG) zero(RG) zero(RG);
    zero(RG) zero(RG) zero(RG) zero(RG) zero(RG) zero(RG);
    zero(RG) zero(RG) zero(RG) zero(RG) zero(RG) zero(RG); 
    zero(RG) zero(RG) zero(RG) zero(RG) zero(RG) zero(RG)
]

@assert act_on_matrix(M, one(Σ), alphabet_permutation_) == M
for σ in Σ
    for τ in Σ
        @assert act_on_matrix(M, σ*τ, alphabet_permutation_) == act_on_matrix(act_on_matrix(M, σ, alphabet_permutation_), τ, alphabet_permutation_)
        # @assert act_on_matrix(M, σ*τ, alphabet_permutation_) == act_on_matrix(act_on_matrix(M, τ, alphabet_permutation_), σ, alphabet_permutation_) # for left action
    end
end
#############################################################################################

function weyl_symmetrize_matrix(
    M::AbstractMatrix{<:AlgebraElement},
    Σ, # the symmetry group (either symmetric group or wreath product)
    op,
    n::Integer
)
    RG = parent(first(M))
    G = parent(first(basis(RG)))
    S = gens(G)

    alphabet_permutation_ = alphabet_permutation(alphabet(G), Σ, op, n)

    result = [zero(RG) for i in eachindex(S), j in eachindex(S)]
    for σ in Σ
        result += act_on_matrix(M, σ, alphabet_permutation_)
    end

    return result
end

Δ₁_emb_symmetrized = let n = 4
    Σ = PermutationGroups.SymmetricGroup(n)
    weyl_symmetrize_matrix(Δ₁_emb, Σ, sln_alphabet_op, n)
end

# The last functions to call:
alpha_beta_adjust(Δ₁⁺, Δ₁⁻, Δ₁_emb_symmetrized)
alpha_beta_adjust(Δ₁⁺, Δ₁⁻, Δ₁_emb_symmetrized, false)

Δ₁⁺[1,1]
Δ₁⁻[1,1]
Δ₁_emb_symmetrized[1,1]

M = 6*Δ₁⁺+6*Δ₁⁻-Δ₁_emb_symmetrized
length(basis(RG_prime))*144

M[1,2]
M[2,2]
M[3,3]

using LinearAlgebra

diff = let M = 12*Δ₁⁺+36*Δ₁⁻-Δ₁_emb_symmetrized
    M_max_I = let
        M_max_I = [0.0*one(RG_prime) for i in 1:20, j in 1:20]
        for i in 1:12
            M_max_I[i,i] = Float64(M[i,i](one(sl4)))*one(RG_prime)
        end
        M_max_I
    end

    M_diff = M-M_max_I

    l1_norm = sum(x -> norm(x, 1), M_diff)
    M_diff[1,1](one(sl4))-l1_norm
end
