# using Pkg
# Pkg.activate(@__DIR__)

# using LinearAlgebra
# ENV["JULIA_NUM_THREADS"] = 4
# LinearAlgebra.BLAS.set_num_threads(2)

include("starAlgebras.jl")

function foxDerivative(relatorWord, generatorId)
    relatorWordLength = length(relatorWord)
    result = []
    multiplier = []
 
    for i in 1:relatorWordLength
        if relatorWord[i] == (2*generatorId-1) # generator
            append!(result, [[copy(multiplier),1]])
        elseif relatorWord[i] == 2*generatorId # inverse of a generator
            tempMultiplier = copy(multiplier)
            append!(tempMultiplier, relatorWord[i])
            append!(result, [[copy(tempMultiplier),-1]])
        end
    
        append!(multiplier, relatorWord[i])
    end
 
    return result
 end
 
 function projectEntryFromFreeGroupRing(freeGroupRingElement, groupGenerators, RG, G)
    result = RG(0)
 
    for i in 1:length(freeGroupRingElement)
        summand = one(G)
        for j in 1:length(freeGroupRingElement[i][1])
            if freeGroupRingElement[i][1][j]%2 == 1
                summand *= groupGenerators[(floor(Int, (freeGroupRingElement[i][1][j]+1)/2))]
            else
                summand *= inv(groupGenerators[(floor(Int, freeGroupRingElement[i][1][j]/2))])
            end
        end
        if freeGroupRingElement[i][2] == 1
            result += RG(summand)
        else
            result -= RG(summand)
        end
    end
 
    return result
 end
 
 function jacobianMatrixEncoded(relations, generators)
    relationsNumber = length(relations)
    generatorsNumber = length(generators)
 
    result = [[[0],1] for i in collect(1:relationsNumber*generatorsNumber)]
    result = reshape(result, relationsNumber, generatorsNumber)
 
    for i in 1:relationsNumber
        for j in 1:generatorsNumber
            relator = relations[i].first*inv(relations[i].second)
            result[i,j] = foxDerivative(relator.word, j)
        end
    end
 
    return result
 end
 
 function jacobianMatrixDecoded(J, generators, RG, G)
    relationsNumber = size(J)[1]
    generatorsNumber = size(J)[2] # which is equal to length(generators) by the way
    result = [RG(0) for i in collect(1:relationsNumber*generatorsNumber)]
    result = reshape(result, relationsNumber, generatorsNumber)
 
    for i in 1:relationsNumber
        for j in 1:generatorsNumber
            result[i,j] = projectEntryFromFreeGroupRing(J[i,j], generators, RG, G)
        end
    end
 
    return result
 end
 
 function suitableGroupRing(G, generators, JEncoded)
    summandElements = [one(G)]
    for i in 1:size(JEncoded)[1]
        for j in 1:size(JEncoded)[2]
            for k in 1:length(JEncoded[i,j])
                summand = one(G)
                for l in 1:length(JEncoded[i,j][k][1])
                    if JEncoded[i,j][k][1][l]%2 == 1
                        summand *= generators[(floor(Int, (JEncoded[i,j][k][1][l]+1)/2))]
                    else
                        summand *= inv(generators[(floor(Int, JEncoded[i,j][k][1][l]/2))])
                    end
                end
                append!(summandElements, [summand])
            end
        end
    end
    support = unique([generators; summandElements])

    return groupRing(G, support)
 end

 function jacobianMatrix(G, JEncoded, generators, RG)
    return jacobianMatrixDecoded(JEncoded, generators, RG, G)
 end
 
 function starOfMatrixOverGroupRing(M)
    result = reshape(copy(M), size(M)[2], size(M)[1])
    for i in 1:size(M)[2]
        for j in 1:size(M)[1]
            result[i,j] = StarAlgebras.star(M[j,i])
        end
    end

    return result
 end

 function D₀(G, generators, RG)
    result = [RG(0) for i in 1:length(generators)]
    result = reshape(result, length(generators), 1)
    for i in 1:length(generators)
        result[i,1] = RG(generators[i])-RG(one(G))
    end

    return result
 end

 function Δ₁⁺(G, JEncoded, generators, RG)
    D₁ = jacobianMatrix(G, JEncoded, generators, RG)

    return starOfMatrixOverGroupRing(D₁)*D₁
 end

 function Δ₁⁻(G, generators, RG)
    D₀x = D₀(G, generators, RG)
    return D₀x*starOfMatrixOverGroupRing(D₀x)
 end

 function Δ₁(G, JEncoded, generators, RG)
    return Δ₁⁺(G, JEncoded, generators, RG)+Δ₁⁻(G, generators, RG)
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
 
 function differentials(G)
    generators = gens(G)
    jacobianMatrixEncodedx = jacobianMatrixEncoded(G.relations, generators)
    RG = suitableGroupRing(G, generators, jacobianMatrixEncodedx)

    D₁ = jacobianMatrix(G, jacobianMatrixEncodedx, generators, RG)

    @info "D₁ = Jacobian"
    printMatrix(D₁)
    # @info D₁
    @info "D₁*"
    printMatrix(starOfMatrixOverGroupRing(D₁))
    # @info "typeof Jacobian:"
    # @info typeof(D₁)

    D₀x = D₀(G, generators, RG)

    @info "D₀"
    printMatrix(D₀x)
    @info "D₀*"
    printMatrix(starOfMatrixOverGroupRing(D₀x))
    
    Δ₁⁺x = Δ₁⁺(G, jacobianMatrixEncodedx, generators, RG)

    @info "Δ₁⁺ = D₁*D₁"
    printMatrix(Δ₁⁺x)

    Δ₁⁻x = Δ₁⁻(G, generators, RG)

    @info "Δ₁⁻ = D₀D₀*"
    printMatrix(Δ₁⁻x)

    Δ₁x = Δ₁(G, jacobianMatrixEncodedx, generators, RG)

    @info "Δ₁ = Δ₁⁺+Δ₁⁻"
    printMatrix(Δ₁x)
 end
