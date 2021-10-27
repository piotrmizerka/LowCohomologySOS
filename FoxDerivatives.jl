using Pkg
Pkg.activate(@__DIR__)

using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = 4
LinearAlgebra.BLAS.set_num_threads(2)

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
            # result[i,j] = projectEntryFromFreeGroupRing(foxDerivative(relator.word, j), generators, RG, G)
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
 
 function jacobianMatrix(G)
    relations = G.relations
    generators = gens(G)

    # maxRelatorWordLength = maximum([length(relation.first.word) for relation in relations])
    # RG = groupRing(G, ceil(Int, maxRelatorWordLength/2)+1)
 
    jacobianEncoded = jacobianMatrixEncoded(relations, generators)

    summandElements = [one(G)]
    # summandElements = [1,1]
    # @info summandElements
    append!(summandElements, [one(G)])
    for i in 1:size(jacobianEncoded)[1]
        for j in 1:size(jacobianEncoded)[2]
            for k in 1:length(jacobianEncoded[i,j])
                summand = one(G)
                for l in 1:length(jacobianEncoded[i,j][k][1])
                    if jacobianEncoded[i,j][k][1][l]%2 == 1
                        summand *= generators[(floor(Int, (jacobianEncoded[i,j][k][1][l]+1)/2))]
                    else
                        summand *= inv(generators[(floor(Int, jacobianEncoded[i,j][k][1][l]/2))])
                    end
                end
                # @info summand
                # @info summandElements
                append!(summandElements, [summand])
            end
        end
    end
    summandElements = unique(summandElements)

    # @info summandElements

    RG = groupRing(G, summandElements)

    result = jacobianMatrixDecoded(jacobianEncoded, generators, RG, G)
 
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
 
 function starOfMatrixOverGroupRing(M)
    result = reshape(copy(M), size(M)[2], size(M)[1])
    for i in 1:size(M)[2]
        for j in 1:size(M)[1]
            result[i,j] = StarAlgebras.star(M[j,i])
        end
    end

    return result
 end

 function Δ₁⁺(G)
    D₁ = jacobianMatrix(G)
    D₁Star = starOfMatrixOverGroupRing(D₁)
 
    # printMatrix(D₁)
 
    result = D₁Star*D₁ # for this moment - just a test of matrix multiplication - works for square matrices only
 
    return result
 end
 
 let
    # A = Alphabet([:x, :X, :y, :Y], [2, 1, 4, 3])
    # F = FreeGroup(A)
    # x,y = Groups.gens(F)
    # ε = one(F)
    # G = FPGroup(F, [x*y => ε] ) # proper matrix: [1 x]
    # G = FPGroup(F, [x*y => ε, y*x => ε] ) # proper matrix: [1 x; X 1]
    # G = FPGroup(F, [x*y*x^(-1)*y^(-1) => ε] ) # proper matrix: [1-y x-1]
    # G = FPGroup(F, [x*y => y*x] ) # proper matrix: as above
 
    A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    x,y,z = Groups.gens(F)
    ε = one(F)
    G = FPGroup(F, [x*y*z => y*x] ) # proper matrix: [1-y x-1 xy]
    # G = FPGroup(F, [x*y*z => y*x, z => ε] ) # proper matrix: [1-y x-1 xy; 0 0 1]
    
    jacobianMatrixx = jacobianMatrix(G)
    printMatrix(jacobianMatrixx)
    # printMatrix(StarAlgebras.star.(jacobianMatrixx))
    printMatrix(starOfMatrixOverGroupRing(jacobianMatrixx))
 
    Δ₁⁺x = Δ₁⁺(G)
 
    # println("")
    printMatrix(Δ₁⁺x)
 end
 
 CₙJacobianMatrix = let n = 10
    Cₙ = cyclicGroup(n)
    jacobianMatrixx = jacobianMatrix(Cₙ)
 
    @info jacobianMatrixx
 end
 
 SL₃ƵShorterPresentation = let maxRules = 1000
    A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    x,y,z = Groups.gens(F)
    ε = one(F);
    G = FPGroup(F, [x^3 => ε, y^3 => ε, z^2 => ε, (x*z)^3 => ε, (y*z)^3 => ε, 
                    (x^(-1)*z*x*y)^2 => ε, (y^(-1)*z*y*x)^2 => ε, (x*y)^6 => ε], maxrules = maxRules )
    
    G
 end
 
 SL₃ƵJacobianShorterPresentation = let
    jacobianMatrixx = jacobianMatrix(SL₃ƵShorterPresentation)
 
    # @info jacobianMatrixx
    printMatrix(jacobianMatrixx)
 end
 
 SL₃ƵElementaryMatrixPresentation = let maxRules = 1000
    A = Alphabet([:e12, :E12, :e21, :E21, :e13, :E13, :e31, :E31, :e23, :E23, :e32, :E32], [2, 1, 4, 3, 6, 5, 8, 7, 10, 9, 12, 11])
    F = FreeGroup(A)
    e12, e21, e13, e31, e23, e32 = Groups.gens(F)
    ε = one(F);
    G = FPGroup(F, [e12*e13 => e13*e12, e12*e32 => e32*e12, e13*e23 => e23*e13, e23*e21 => e21*e23, e21*e31 => e31*e21, e31*e32 => e32*e31,
                    e12*e23 => e13*e23*e12, e13*e32 => e12*e32*e13, e21*e13 => e23*e13*e21, e23*e31 => e21*e31*e23, 
                    e31*e12 => e32*e12*e31, e32*e21 => e31*e21*e32,
                    (e12*e21^(-1)*e12)^4 => ε], maxrules = maxRules )
    
    G
 end
 
 SL₃ƵJacobianElementaryMatrixPresentation = let
    jacobianMatrixx = jacobianMatrix(SL₃ƵElementaryMatrixPresentation)
 
    # @info jacobianMatrixx
    printMatrix(jacobianMatrixx)
 end