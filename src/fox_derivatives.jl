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

    F = parent(u)
    xᵢ = Groups.gens(F, i)
    current_multiplier = one(F)
    result = zero(RF)
    for j in 1:length(word(u))
        if F(word(u)[j:j]) == xᵢ
            result += RF(current_multiplier)
        elseif F(word(u)[j:j]) == inv(xᵢ)
            result -= RF(current_multiplier*inv(xᵢ))
        end
        current_multiplier *= F(word(u)[j:j])
    end

    return result
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
    sizehint!(half_basis, length(relations) * maximum(length ∘ word, relations))

    for rel in relations
        for k in 1:length(word(rel))
            for l in k:length(word(rel))
                append!(half_basis, [F(word(rel)[k:l])])
            end
        end
    end

    return group_ring(F, half_basis, additive_only = true)
end
