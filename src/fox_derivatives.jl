function d₀(RG, generators)
    result = [RG(g)-one(RG) for g in generators]
    return reshape(result, length(generators), 1)
end

# h is intended to be a homomorphism from a free group to G
function embed_to_group_ring(X::AlgebraElement, RG::StarAlgebra, h::Function)
    length(supp(X)) == 0 && return zero(RG)
    return sum(X(g)*RG(h(g)) for g in supp(X))
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
        d = div(length(word(u)),2)
        p = parent(u)(word(u)[1:d])
        s = parent(u)(word(u)[d+1:end])

        return fox_derivative(RF, p, i) + RF(p)*fox_derivative(RF, s, i)
    end
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
    F = parent(rand(relations))
    half_basis = [one(F)]

    function relation_append_basis(u::FPGroupElement)
        for k in 1:length(word(u))
            for l in k:length(word(u))
                append!(half_basis, [F(word(u)[k:l])])
            end
        end
    end

    for r in relations
        relation_append_basis(r)
    end

    return group_ring(F, half_basis)
end
