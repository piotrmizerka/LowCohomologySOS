# We want to find alpha and beta such that
# αd₁^*d₁+βd₀d₀^* = symmetrized induced Laplacians from SL(3,Z) to SL(4,Z)

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
    i::Groups.Homomorphism;
    half_radius::Integer # stands for the half radius for the group ring of G' being the support of the Laplacian for G' (to be computed separately)
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

    RG_prime = LowCohomologySOS.group_ring(G_prime, half_radius)
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

    SLₙℤ = SL(n, UInt8)
    SLₘℤ = SL(m, UInt8)

    _idx(k) = ((i,j) for i in 1:k for j in 1:k if i≠j)

    inds_S_SLₙℤ = Dict( 
        let eij = MatrixGroups.ElementaryMatrix{n}(i,j, UInt8(1))
            SLₙℤ([alphabet(SLₙℤ)[eij]])
        end => (i,j)
        for (i,j) in _idx(n)
    )
    S_SLₘℤ = Dict((i,j) =>
        let eij = MatrixGroups.ElementaryMatrix{m}(i,j, UInt8(1))
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
M_emb = embed_matrix(M, i, half_radius = 2)
#############################################################################################

function laplacian_embedding(n::Integer, m::Integer)
    # TODO ??
end

# Laplacian embedding from SL(3,Z) to SL(4,Z) ###############################################
const half_radius = 2

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

Δ₁, Iₙ = LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)

Δ₁_emb = embed_matrix(Δ₁, i, half_radius = 2)
#############################################################################################