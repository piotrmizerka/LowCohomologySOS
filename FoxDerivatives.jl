# using Oscar

include("starAlgebras.jl")

function foxDerivative(RF::StarAlgebra, u::FPGroupElement, i::Integer)
    isone(u) && return RF(0)

    if length(word(u)) == 1
        if u == Groups.gens(parent(u), i)
            return one(RF)
        end
        if u == inv(Groups.gens(parent(u), i))
            return -RF(u)
        else
            return RF(0)
        end
    else
        half = ceil(Int,length(word(u))/2)

        # @info half

        p = parent(u)(word(u)[1:half])
        s = parent(u)(word(u)[(half+1):end])

        # @info p
        # @info s

        return foxDerivative(RF, p, i) + RF(p)*foxDerivative(RF, s, i)
    end
end

# h is intended to be a homomorphism from a free group to G
function embedToGroupRing(h::Function, RG::StarAlgebra, x::AlgebraElement)
    if length(supp(x)) == 0
        return RG(0)
    end

    return sum(x(g)*RG(h(g)) for g in supp(x))
end

# For Fox derivatives in the free group ring. NOT WORKING AS INTENDED!
# function suitableGroupRing(relations)
#     F = parent(relations[1])
#     summandElements = [one(F), one(F)]
#     for r in relations
#         for i in 2:length(r.word)
#             # append!(summandElements, [F(r.word[1:i])])
#             append!(summandElements, [F(r.word[(i-1):length(r.word)])])
#         end
#     end
#     halfBasis = unique([Groups.gens(F); inv.(Groups.gens(F)); summandElements])
    
#     return groupRing(F, halfBasis)
# end

# function suitableGroupRing(J::Matrix{<:StarAlgebra}, h::Function) PROBLEMS WITH TYPES
function suitableGroupRing(J, h::Function)
    RF = parent(J[1,1])
    G = parent(h(RF.basis[1]))
    summandElements = [one(G)]
    for Jᵢⱼ in J
        for x in supp(Jᵢⱼ)
            # @info typeof(summandElements[1])
            # @info typeof(h(x))

            append!(summandElements, [h(x)])
        end
    end
    halfBasis = unique(summandElements)
    
    return groupRing(G, halfBasis)
end

# Jacobian matrix in the free group ring
# Relations is an array of relators which are elements of the free group
function jacobianMatrix(relations)
    F = parent(relations[1])

    # @info F

    relationsNumber = length(relations)
    generatorsNumber = length(Groups.gens(F))

    # RF = suitableGroupRing(relations)
    halfRadius = ceil(Int, maximum([length(r.word) for r in relations])/2)

    # @info halfRadius

    RF, halfBasis  = groupRing(F, halfRadius)

    # @info length(RF.basis)
 
    result = [RF(0) for i in 1:relationsNumber*generatorsNumber]
    result = reshape(result, relationsNumber, generatorsNumber)
 
    for i in 1:relationsNumber
        for j in 1:generatorsNumber
            # @info typeof(result[i,j])
            # @info typeof(foxDerivative(RF, relations[i], j))

            result[i,j] = foxDerivative(RF, relations[i], j)
        end
    end
 
    return result
end

# Jacobian matrix in the proper ring. 
# Requires group homomorphism h from the free group to the considered group
# and precomputed free group Jacobian matrix J of the considered group.
function jacobianMatrix(J, h::Function, RG)
    relationsNumber = size(J)[1]
    generatorsNumber = size(J)[2]
    
    # Determine the Jacobian in RG
    result = [RG(0) for i in 1:relationsNumber*generatorsNumber]
    result = reshape(result, relationsNumber, generatorsNumber)
 
    for i in 1:relationsNumber
        for j in 1:generatorsNumber
            result[i,j] = embedToGroupRing(h, RG, J[i,j])
        end
    end
 
    return result
end

function starOfMatrixOverGroupRing(M)
    RG = parent(M[1,1])
    result = reshape(copy(M), size(M)[2], size(M)[1])
    for i in 1:size(M)[2]
        for j in 1:size(M)[1]
            # result[i,j] = StarAlgebras.star(M[j,i])
            # @info M[i,j].coeffs[length(RG.basis)]*RG(inv(RG.basis[1]))

            result[i,j] = sum(M[j,i].coeffs[k]*RG(inv(RG.basis[k])) for k in 1:length(RG.basis))
            # result[i,j] = M[i,j].coeffs[1]*RG(inv(RG.basis[1]))
            # result[i,j] = RG(0)

            # @info "i,j:"
            # @info i 
            # @info j

            # for k in 1:length(RG.basis)
                # @info k
                # @info M[i,j].coeffs[k]
                # @info RG(inv(RG.basis[k]))

                # result[i,j] += M[i,j].coeffs[k]*RG(inv(RG.basis[k]))
            # end
        end
    end

    return result
end

function D₀(G, RG, generators)
    result = [RG(0) for i in 1:length(generators)]
    result = reshape(result, length(generators), 1)
    for i in 1:length(generators)
        result[i,1] = RG(generators[i])-RG(one(G))
    end

    return result
end

function printMatrix(M)
    for i in 1:size(M, 1)
        for j in 1:size(M, 2)
            print(M[i,j])
            print("   ")
        end
        println("")
    end
    println("")
end
