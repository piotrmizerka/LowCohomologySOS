using StarAlgebras
using AbstractAlgebra
using Groups
using Test

Base.:/(G::Groups.AbstractFPGroup, rels::AbstractVector{<:FPGroupElement}) = FPGroup(G, [r=>one(G) for r in rels])
StarAlgebras.star(g::Groups.GroupElement) = inv(g)
∗ = StarAlgebras.star


# Group ring with basis the whole ball
function groupRing(G, supportSize::Int64, starMultiplication = false)
    S = Groups.gens(G)
    S = unique([S; inv.(S)])
    ID = one(G)
    Ball, sizes = Groups.wlmetric_ball(S, ID, radius = 2*supportSize)

    # @info "Ball size:"
    # @info length(Ball)

    b = StarAlgebras.Basis{UInt32}(Ball)
    tmstr = StarAlgebras.CachedMTable{starMultiplication}(b, table_size = (sizes[supportSize], sizes[supportSize]))
    RG = StarAlgebra(G, b, tmstr)

    return RG
end

# Group ring with basis given by prescribed support
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
    groupRingBasis = StarAlgebras.Basis{UInt32}(basisElements)
    tmstr = StarAlgebras.CachedMTable{starMultiplication}(groupRingBasis, table_size = (length(inverseClosedSupport), length(inverseClosedSupport)))
    RG = StarAlgebra(G, groupRingBasis, tmstr)

    return RG
end

# Saves a matrix MRGOld over the group ring MRGOld as a matrix over the group ring XRGNew (both rings are for G)
function changeUnderlyingGroupRing(MRGOld, RGOld, RGNew, G)
    result = [RGNew(0) for i in 1:length(MRGOld)]
    result = reshape(result, size(MRGOld)[1], size(MRGOld)[2])
    for i in 1:size(MRGOld)[1]
        for j in 1:size(MRGOld)[2]
            for k in 1:length(MRGOld[i,j].coeffs)
                if MRGOld[i,j].coeffs[k] != 0
                    result[i,j] += MRGOld[i,j].coeffs[k]*RGNew(RGOld.basis[k])
                end
            end
        end
    end

    return result
end

function cyclicGroup(n)
    A = Alphabet([:a, :A, :b, :B], [2,1,4,3])
    F = FreeGroup(A)
    a, b = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e, b => e])

    return Cₙ
end
