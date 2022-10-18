function d₀(RG, generators)
    result = [RG(g) - one(RG) for g in generators]
    return reshape(result, length(generators), 1)
end

# h is intended to be a homomorphism from a free group to G
function embed(h, X::AlgebraElement, RG::StarAlgebra)
    S = supp(X)
    length(S) == 0 && return zero(RG, eltype(X))
    return sum(X(s) * RG(h(s)) for s in S)
end

function fox_derivative(RF::StarAlgebra, u::FPGroupElement, i::Integer)
    @assert parent(u) === StarAlgebras.object(RF)

    isone(u) && return zero(RF)

    if length(word(u)) == 1
        g = Groups.gens(parent(u), i)
        if u == g
            return one(RF)
        elseif u == inv(g)
            return -RF(inv(g))
        else
            return zero(RF)
        end
    else
        d = div(length(word(u)), 2)
        p = parent(u)(word(u)[1:d])
        s = parent(u)(word(u)[d+1:end])

        return fox_derivative(RF, p, i) + RF(p) * fox_derivative(RF, s, i)
    end
end

function fox_derivative(F::FreeGroup, u::FPGroupElement, i::Integer)
    coeffs = Int[]
    sizehint!(coeffs, length(word(u)))
    elts = eltype(F)[]
    sizehint!(elts, length(word(u)))

    @assert alphabet(F) == alphabet(parent(u))
    A = alphabet(F)
    li = let g = Groups.gens(F, i)
        @assert length(word(g)) != 0
        first(word(g))
    end
    current = one(F)

    for letter in word(u)
        if letter == li
            push!(coeffs, 1)
            push!(elts, F(copy(word(current))))
        elseif inv(A, letter) == li
            push!(coeffs, -1)
            push!(elts, current*inv(F([li])))
        end
        push!(word(current), letter)
    end
    return coeffs, elts
end

# Jacobian matrix in the free group ring.
# Relations is an array of relators which are elements of the free group
function jacobian_matrix(relations)
    @assert !isempty(relations)
    F = parent(first(relations))
    RF = suitable_group_ring(relations)

    jac = [fox_derivative(RF, r, j) for r in relations, j in 1:Groups.ngens(F)]

    return jac
end

# relations is intended to contain elements of a free group
function suitable_group_ring(relations)
    @assert !isempty(relations)
    F = parent(first(relations))

    half_basis = [one(F)]
    for s in gens(F)
        append!(half_basis, [s])
    end
    sizehint!(half_basis, length(relations) * maximum(length ∘ word, relations))

    for rel in relations
        for k in 1:length(word(rel))
            for l in k:length(word(rel))
                append!(half_basis, [F(word(rel)[k:l])])
            end
        end
    end

    half_basis = unique!(half_basis)

    return group_ring(F, half_basis, additive_only = true)
end
