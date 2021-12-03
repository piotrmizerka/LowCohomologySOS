Base.:/(G::Groups.AbstractFPGroup, rels::AbstractVector{<:FPGroupElement}) = FPGroup(G, [r=>one(G) for r in rels])
StarAlgebras.star(g::Groups.GroupElement) = inv(g)
StarAlgebras.star(A::AbstractAlgebra.Generic.MatAlgElem{Int64}) = inv(A)
∗ = StarAlgebras.star

# UWAGI OGOLNE:
# RG.basis zamienic na basis(RG)
# Revise - sledzi zmiany - po uaktualnieniach
# @edit zero(RG) - tak mozna zobaczyc jak funkcja jest zdefiniowa

# Group ring with basis the whole ball; returns halfbasis as well
function groupRing(G, halfRadius::Integer, starMultiplication = false)
    S = Groups.gens(G)
    S = unique([S; inv.(S)])
    ID = one(G)
    Ball, sizes = Groups.wlmetric_ball(S, ID, radius = 2*halfRadius)
    b = StarAlgebras.Basis{UInt32}(Ball)
    tmstr = StarAlgebras.CachedMTable{starMultiplication}(b, table_size = (sizes[halfRadius], sizes[halfRadius]))
    RG = StarAlgebra(G, b, tmstr)

    # Nie trzeba zrwaca kuli - ta informacja jest w RG
    return RG, Ball[1:sizes[halfRadius]]
end

# Group ring with basis given by prescribed support
function groupRing(G, halfBasis, starMultiplication = false)
    # inverseClosed - lepiej inverse_closed - uwaga styl
    starClosedSupport = unique!([halfBasis; star.(halfBasis)]) # ! modyfikuje in place (zysk jednej alokacji)
    basisElements = let basis = starClosedSupport, f = starMultiplication ? star : identity
        unique!([basis;[f(a)*b for a in basis for b in basis]])
    end
    groupRingBasis = StarAlgebras.Basis{UInt32}(basisElements)
    tmstr = StarAlgebras.CachedMTable{starMultiplication}(groupRingBasis, table_size = (length(starClosedSupport), length(starClosedSupport)))
    RG = StarAlgebra(G, groupRingBasis, tmstr)

    return RG
end

# Saves a matrix MRGOld over the group ring MRGOld as a matrix over the group ring RGNew (both rings are for G)
function changeUnderlyingGroupRing(MRGOld, RGOld, RGNew)
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

function changeUnderlyingGroupRing(X::AlgebraElement, RGNew::StarAlgebra)
    @assert StarAlgebras.object(parent(X)) === StarAlgebras.object(RGNew)
    result = RGNew(zero(eltype(X))) #zero(RGNew)
    for x in supp(X)
        # result += X(x)*RGNew(x) may be slow due to allocation
        result[x] += X(x)
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

