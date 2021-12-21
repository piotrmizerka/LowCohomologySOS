function d₀(RG, generators)
    result = reshape([zero(RG) for i in 1:length(generators)], length(generators), 1)
    
    for i in 1:length(generators)
        result[i,1] = RG(generators[i])-one(RG)
    end

    return result
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
    F = parent(rand(relations))
    RF = suitable_group_ring(relations)
    relations_number = length(relations)
    generators_number = length(Groups.gens(F))
 
    result = reshape([zero(RF) for i in 1:relations_number*generators_number], relations_number, generators_number)
 
    for i in 1:relations_number
        for j in 1:generators_number
            result[i,j] = fox_derivative(RF, relations[i], j)
        end
    end
 
    return result
end

function print_matrix(M)
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            print(M[i,j])
            print("   ")
        end
        println("")
    end
    println("")
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
