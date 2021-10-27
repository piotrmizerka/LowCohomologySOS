using StarAlgebras
using AbstractAlgebra
using Groups
using Test

Base.:/(G::Groups.AbstractFPGroup, rels::AbstractVector{<:FPGroupElement}) = FPGroup(G, [r=>one(G) for r in rels])
StarAlgebras.star(g::Groups.GroupElement) = inv(g)
∗ = StarAlgebras.star

# Symmetric group ring
function freeGroupRing(generatorsNumber, radius)
    println("\n*********FREE GROUP RING************\n")

    G = FreeGroup(generatorsNumber)
    S = gens(G)
    S = unique([S; inv.(S)])
    Bᵣ, sizes = Groups.wlmetric_ball(S, one(G), radius = 2 * radius)
    b = StarAlgebras.Basis{UInt32}(Bᵣ)
    tmstr = StarAlgebras.MTable{true}(b, table_size = (sizes[radius], sizes[radius]))
    RG = StarAlgebra(G, b, tmstr)

    println("\nGROUP DEFINING THE MODULE:\n",RG.object)
    println("\nBASIS:\n",RG.basis)
    println("\nMULTIPLICATION TABLE:\n",RG.mstructure)
end;

# Symmetric group ring - not working
function symmetricGroupRing(n, displayMode = false)
    G = SymmetricGroup(n)
    b = StarAlgebras.Basis{UInt32}(collect(G)) # this causes problems in the definition of tmstr

    @info b

    tmstr = StarAlgebras.MTable{true}(b, table_size = (factorial(n)/2, factorial(n)/2))
    RG = StarAlgebra(G, b, tmstr)

    if displayMode
        println("\n*********SYMMETRIC GROUP RING************\n")
        println("GROUP DEFINING THE MODULE:\n",RG.object)
        println("\nBASIS:\n",RG.basis)
        println("\nMULTIPLICATION TABLE:\n",RG.mstructure)
    else
        return RG
    end
end;

# Group ring of G_1 form the article of Fujiwara from 2017
function G1GroupRing(halfBasisLength = 1, displayMode = false)
    A = Alphabet([:a, :A, :b, :B, :c, :C], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    a,b,c = Groups.gens(F)
    ε = one(F);
    G1 = FPGroup(F, [a^3 => ε, b^3 => ε, c^3 => ε, 
                     (a*b)^2 => b*a, (b*c)^2 => c*b, (c*a)^2 => a*c ], maxrules = 238)
    S = Groups.gens(G1)
    S = unique([S; inv.(S)])
    ID = one(G1)
    Bᵣ, sizes = Groups.wlmetric_ball(S, ID, radius = 2 * halfBasisLength)
    b = StarAlgebras.Basis{UInt32}(Bᵣ)
    tmstr = StarAlgebras.MTable{true}(b, table_size = (sizes[halfBasisLength], sizes[halfBasisLength]))
    RG1 = StarAlgebra(G1, b, tmstr)

    if displayMode
        println("GROUP:\n", G1)
        println("\nBALL OF RADIUS ", 2*halfBasisLength, ":\n", Bᵣ)
        println("\nSIZES OF BALLS FOR SUBSEQUENT RADII:\n",sizes)
        println("\nGROUP RING OF G1:\n",RG1)
        println("\nHAS BASIS? ", isdefined(RG1,:basis))
        println("\nZERO OF THE STAR ALGEBRA: ",zero(RG1))
        println("\nONE OF THE STAR ALGEBRA: ",one(RG1))
        println("\nONE OF THE STAR ALGEBRA (USING DIFFERENT CALL): ",RG1(1))
    else
        return RG1, ID
    end
end;

# Group ring with basis the whole ball
function groupRing(G, supportSize::Int64, starMultiplication = false)
    S = Groups.gens(G)
    S = unique([S; inv.(S)])
    ID = one(G)
    Ball, sizes = Groups.wlmetric_ball(S, ID, radius = 2*supportSize)

    # @info sizes

    b = StarAlgebras.Basis{UInt32}(Ball)
    tmstr = StarAlgebras.CachedMTable{starMultiplication}(b, table_size = (sizes[supportSize], sizes[supportSize]))
    RG = StarAlgebra(G, b, tmstr)
    # RG = StarAlgebra(G, b)

    return RG
end;

# Group ring with basis given by prescribed support - TODO
function groupRing(G, support::AbstractVector{<:FPGroupElement}, starMultiplication = false)
    basisElementsToAppend = []
    inverseClosedSupport = unique([support; inv.(support)])
    for i in 1:length(inverseClosedSupport)
        for j in 1:length(inverseClosedSupport)
            g = inverseClosedSupport[i]*inverseClosedSupport[j]
            if !(g in inverseClosedSupport)
                append!(basisElementsToAppend, [g])
            end
        end
    end
    basisElementsToAppend = unique(basisElementsToAppend)
    basisElements = copy(inverseClosedSupport)
    for g in basisElementsToAppend
        append!(basisElements, [g])
    end

    # @info basisElements

    groupRingBasis = StarAlgebras.Basis{UInt32}(basisElements)
    tmstr = StarAlgebras.CachedMTable{starMultiplication}(groupRingBasis, table_size = (length(inverseClosedSupport), length(inverseClosedSupport)))
    RG = StarAlgebra(G, groupRingBasis, tmstr)

    return RG
end

function cyclicGroup(n)
    A = Alphabet([:a, :A, :b, :B], [2,1,4,3])
    F = FreeGroup(A)
    a, b = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e, b => e])

    return Cₙ
end;

# Cyclic group ring
function cyclicGroupRing(n)
    Cₙ = cyclicGroup(n)
    ID = one(Cₙ)
    RCₙ = groupRing(Cₙ, n)
    return RCₙ, ID
end;
