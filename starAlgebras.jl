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

# Symmetric group ring
function symmetricGroupRing(n, displayMode = false)
    G = SymmetricGroup(n)
    RG = StarAlgebra(G, StarAlgebras.Basis{UInt8}(collect(G)))

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

# Cyclic group ring
function cyclicGroupRing(n)
    A = Alphabet([:a, :A, :b, :B], [2,1,4,3])
    F = FreeGroup(A)
    a, b = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e, b => e])
    S = Groups.gens(Cₙ)
    S = unique([S; inv.(S)])
    ID = one(Cₙ)
    Bᵣ, sizes = Groups.wlmetric_ball(S, ID, radius = n)
    b = StarAlgebras.Basis{UInt32}(Bᵣ)
    basis = StarAlgebras.Basis{UInt8}(b)
    RCₙ = StarAlgebra(Cₙ, basis)

    return RCₙ, ID
end;
