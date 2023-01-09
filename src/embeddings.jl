# Suppose we have a matrix M indexed by the generators S of a group G which is embedded in G' such that S
# is a subset of S', the generating set of G' (the embedding is i:G-->G', i(S)⊆S'). 
# Then the embedding below takes M to M', indexed by S', as follows: M'ₛₜ = Mₛₜ for s,t∈S, and M'ₛₜ = 0 otherwise.
function embed_matrix(
    M::AbstractMatrix{<:AlgebraElement},
    i::Groups.Homomorphism, # i is intended to be an embedding
    RG_prime::StarAlgebra # we must provide the same underlying group ring (it has to be the group ring of i.target)
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

    result = [zero(RG_prime) for i in eachindex(S_prime), j in eachindex(S_prime)]

    for s in S
        for t in S
            result[S_prime_idies[i(s)],S_prime_idies[i(t)]] = embed(i, M[S_idies[s],S_idies[t]], RG_prime)
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