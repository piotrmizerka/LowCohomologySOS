using StarAlgebras
using AbstractAlgebra
using Groups
using Test
# using PropertyT

Base.:/(G::Groups.AbstractFPGroup, rels::AbstractVector{<:FPGroupElement}) = FPGroup(G, [r=>one(G) for r in rels])
StarAlgebras.star(g::Groups.GroupElement) = inv(g)
∗ = StarAlgebras.star

# Symmetric group ring
function S3GroupRing()
    G = SymmetricGroup(3)

    # println(G)

    b = StarAlgebras.Basis{UInt8}(collect(G))

    # println(b)

    RG = StarAlgebra(G, b)

    # println(RG)
    # println("\nHAS BASIS? ", isdefined(RG,:basis))
    # println("\nCOERCING SCALARS:\n",RG(-9.1))

    return RG
end

# Group ring of G_1 form the article of Fujiwara from 2017
function G1GroupRing()
    A = Alphabet([:a, :A, :b, :B, :c, :C], [2, 1, 4, 3, 6, 5])
    F = FreeGroup(A)
    a,b,c = Groups.gens(F)
    ε = one(F);
    G1 = FPGroup(F, [a^3 => ε, b^3 => ε, c^3 => ε, 
                     (a*b)^2 => b*a, (b*c)^2 => c*b, (c*a)^2 => a*c ], maxrules = 238)
    # println("GROUP:\n", G1)

    S = Groups.gens(G1)
    S = unique([S; inv.(S)])
    ID = one(G1)
    RADIUS = 1
    Bᵣ, sizes = Groups.wlmetric_ball(S, ID, radius = 2 * RADIUS)
    # println("\nBALL OF RADIUS ", 2*RADIUS, ":\n", Bᵣ)
    # println("\nSIZES OF BALLS FOR SUBSEQUENT RADII:\n",sizes)

    b = StarAlgebras.Basis{UInt32}(Bᵣ)
    tmstr = StarAlgebras.MTable{true}(b, table_size = (sizes[RADIUS], sizes[RADIUS]))
    RG1 = StarAlgebra(G1, b, tmstr)
    # println("\nGROUP RING OF G1:\n",RG1)
    # println("\nHAS BASIS? ", isdefined(RG1,:basis))
    # println("\nZERO OF THE STAR ALGEBRA: ",zero(RG1))
    # println("\nONE OF THE STAR ALGEBRA: ",one(RG1))
    # println("\nONE OF THE STAR ALGEBRA (USING DIFFERENT CALL): ",RG1(1))

    return RG1, RADIUS, G1
end

# Basis group ring operations - intentionally computation of SOS for G1
function groupRingOperations(RG::StarAlgebra, summandFactorRadius, G::FPGroup)
    S = Groups.gens(G)
    S = unique([S; inv.(S)])
    ID = one(G1)
    Bᵣ, sizes = Groups.wlmetric_ball(S, ID, radius = summandFactorRadius)
    a, b, c = Bᵣ[2:4]
    G = (one(RG) - RG(a))
    H = (one(RG) - RG(b))
    # K = *(RG(c)) * not working...
    K = StarAlgebras.star(RG(c)+RG(c^2))

    println(G^2==2*one(RG)-RG(a)-RG(a^3))
    println(H)
    println(K)
end

# RG1, RADIUS, G1 = G1GroupRing()
# println(RG1.mstructure)
# groupRingOperations(RG1,RADIUS,G1)
# RS3 = symmetricGroup()
# println(RS3.mstructure)