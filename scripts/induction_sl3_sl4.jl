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

# Suppose we have a matrix M indexed by the generators S of a group G which is embedded in G' such that S
# is a subset of S', the generating set of G' (the embedding is i:G-->G', i(S)⊆S'). 
# Then the embedding below takes M to M', indexed by S', as follows: M'ₛₜ = Mₛₜ for s,t∈S, and M'ₛₜ = 0 otherwise.
function embed_matrix(
    M::AbstractMatrix{<:AlgebraElement},
    i::Groups.Homomorphism,
    RG_prime::StarAlgebra # we must provide the same underlying group rin
)
    G = i.source
    G_prime = i.target
    S = gens(G)
    S_prime = gens(G_prime)

    @assert size(M) == (length(S), length(S))

    S_idies = Dict(S[i] => i for i in eachindex(S))
    S_prime_idies = Dict(S_prime[i] => i for i in eachindex(S_prime))

    RG = parent(first(M))
    @assert all(x -> parent(x) === RG, M)
    @assert G == parent(first(basis(RG)))

    # RG_prime = LowCohomologySOS.group_ring(G_prime, half_radius)
    result = [i ≠ j ? zero(RG_prime) : one(RG_prime) for i in eachindex(S_prime), j in eachindex(S_prime)]

    for s in S
        for t in S
            result[S_prime_idies[i(s)],S_prime_idies[i(t)]] = LowCohomologySOS.embed(i, M[S_idies[s],S_idies[t]], RG_prime)
        end
    end

    return result
end

SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)

function sln_slm_embedding(n::Integer, m::Integer)
    @assert n <= m

    SLₙℤ = SL(n, Int8)
    SLₘℤ = SL(m, Int8)

    _idx(k) = ((i,j) for i in 1:k for j in 1:k if i≠j)

    inds_S_SLₙℤ = Dict( 
        let eij = MatrixGroups.ElementaryMatrix{n}(i,j, Int8(1))
            SLₙℤ([alphabet(SLₙℤ)[eij]])
        end => (i,j)
        for (i,j) in _idx(n)
    )
    S_SLₘℤ = Dict((i,j) =>
        let eij = MatrixGroups.ElementaryMatrix{m}(i,j, Int8(1))
            SLₘℤ([alphabet(SLₘℤ)[eij]])
        end
        for (i,j) in _idx(n)
    )
    
    function f(letter_id, SLₙℤ, G)
        if letter_id <= length(gens(SLₙℤ))
            return word(S_SLₘℤ[inds_S_SLₙℤ[SLₙℤ([letter_id])]])
        else
            return word(inv(S_SLₘℤ[inds_S_SLₙℤ[inv(SLₙℤ([letter_id]))]]))
        end
    end

    result = let source = SLₙℤ, target = SLₘℤ
        Groups.Homomorphism(f, source, target, check = false)
    end

    return result
end

# Code below intended for tests ###############################################################
i = sln_slm_embedding(3,4)
sl3 = i.source
sl4 = i.target

e12, e13, e21, e23, e31, e32 = gens(sl3)
i(e12)
i(e13)
i(e21)
i(e23)
i(e31)
i(e32)
i(e12^(-1))
i(e13^(-1))
i(e21^(-1))
i(e23^(-1))
i(e31^(-1))
i(e32^(-1))

RG  = LowCohomologySOS.group_ring(sl3, 2)
M = [
    one(RG) one(RG) zero(RG) zero(RG) RG(gens(sl3,1)) zero(RG);
    RG(gens(sl3,2)) one(RG) one(RG) zero(RG) zero(RG) zero(RG);
    one(RG) one(RG) RG(gens(sl3,3)) zero(RG) zero(RG) zero(RG);
    one(RG) one(RG) zero(RG) zero(RG) RG(gens(sl3,1)) zero(RG);
    one(RG) one(RG) zero(RG) zero(RG) RG(gens(sl3,1)) zero(RG);
    one(RG) one(RG) zero(RG) zero(RG) RG(gens(sl3,1)) zero(RG)
]
RG_prime = LowCohomologySOS.group_ring(sl4,2)
M_emb = embed_matrix(M, i, RG_prime)
#############################################################################################

i = sln_slm_embedding(3,4)
sl4 = i.target

function determine_letter(g)
    @assert length(word(g)) == 1
    
    A = alphabet(parent(g))

    return A[first(word(g))]
end

Δ₁⁺, Δ₁⁻ = let half_radius = 2
    S = gens(sl4)
    S_inv = let s = S
        [s; inv.(s)]
    end
    half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = half_radius)

    F_sl_4_z = FreeGroup(alphabet(sl4))

    quotient_hom = let source = F_sl_4_z, target = sl4
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    N = 4

    elmatrix_gen_dict = Dict(determine_letter(S[i]) => gens(F_sl_4_z, i) for i in eachindex(S))
    e(i,j) = elmatrix_gen_dict[MatrixGroups.ElementaryMatrix{N}(i,j,Int8(1))]

    # Interesting - I didn't know that "==" can return "true" for two elts of different type in Julia:  #######
    # @info typeof(determine_letter(S[1]))
    # @info typeof(MatrixGroups.ElementaryMatrix{N}(1,2))
    # @info typeof(determine_letter(S[1])) == typeof(MatrixGroups.ElementaryMatrix{N}(1,2))
    # @info determine_letter(S[1]) == MatrixGroups.ElementaryMatrix{N}(1,2)
    # @info elmatrix_gen_dict[determine_letter(S[1])]
    #########################################################################################################
    range_as_list = [i for i in 1:N]
    quadruples_total = [(i,j,k,m) for k ∈ 1:N
                            for m ∈ deleteat!(copy(range_as_list), findall(m->m==k,copy(range_as_list)))
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i==m,copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    quadruples_wrong_1 = [(i,j,i,j) for i ∈ 1:N
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    quadruples_wrong_2_inds = []
    for ind in eachindex(quadruples_total)
        (i,j,k,m) = quadruples_total[ind]
        if (i,j) > (k,m)
            append!(quadruples_wrong_2_inds, ind)
        end
    end
    quadruples_wrong_2 = [quadruples_total[ind] for ind in quadruples_wrong_2_inds]
    quadruples = setdiff(setdiff(quadruples_total, quadruples_wrong_1), quadruples_wrong_2)
    triples = [(i,j,k) for i ∈ 1:N
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    
    # The presentation taken from the article of Conder et. al.: https://www.jstor.org/stable/2159559#metadata_info_tab_contents 
    relations = vcat(
        [e(i,j)*e(k,m)*e(i,j)^(-1)*e(k,m)^(-1) for (i,j,k,m) ∈ quadruples],
        [e(i,j)*e(j,k)*e(i,j)^(-1)*e(j,k)^(-1)*e(i,k)^(-1) for (i,j,k) ∈ triples]
    )

    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)

    Δ₁⁺, Δ₁⁻
end


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
        e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1),
    ]

    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)

    Δ₁_emb = embed_matrix(Δ₁, i, RG_prime)
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
    Σ = PermutationGroups.SymmetricGroup(4)
    weyl_symmetrize_matrix(Δ₁_emb, Σ, sln_alphabet_op, n)
end

# The last functions to call:
alpha_beta_adjust(Δ₁⁺, Δ₁⁻, Δ₁_emb_symmetrized)
alpha_beta_adjust(Δ₁⁺, Δ₁⁻, Δ₁_emb_symmetrized, false)
