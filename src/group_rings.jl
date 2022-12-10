Base.:/(G::Groups.AbstractFPGroup, rels::AbstractVector{<:FPGroupElement}) =
    FPGroup(G, [r => one(G) for r in rels])
Base.adjoint(X::AlgebraElement) = StarAlgebras.star(X)
Base.copy(X::AlgebraElement) =
    AlgebraElement(copy(StarAlgebras.coeffs(X)), parent(X))

StarAlgebras.star(g::Groups.AbstractFPGroupElement) = inv(g)
function StarAlgebras.Basis{I}(basis::AbstractVector) where {I}
    # length(basis) <= typemax(I) ||
    #     throw("index type $I is to small for basis of length $(length(basis))")
    # @assert !(eltype(basis) <: Integer)
    # return StarAlgebras.Basis(basis, Dict(b => I(idx) for (idx, b) in pairs(basis)))
    return StarAlgebras.Basis(basis, Dict(b => UInt32(idx) for (idx, b) in pairs(basis)))
end

# Group ring with the basis given by the prescribed support
function group_ring(
    G,
    half_basis;
    star_multiplication = false,
    star_closed=false,
    additive_only = false
)
    star_closed_support = if !star_closed
        unique!([half_basis; star.(half_basis)])
    else
        half_basis
    end

    if !additive_only
        basis_elements =
            let basis = star_closed_support,
                f = star_multiplication ? star : identity

                unique!([basis; [f(a) * b for a in basis for b in basis]])
            end
        group_ring_basis = StarAlgebras.Basis{Int}(basis_elements)
        tmstr = StarAlgebras.CachedMTable{star_multiplication}(
            group_ring_basis,
            table_size = (length(star_closed_support), length(star_closed_support)),
        )
        RG = StarAlgebra(G, group_ring_basis, tmstr)
    else
        group_ring_basis = StarAlgebras.Basis{Int}(star_closed_support)
        RG = StarAlgebra(G, group_ring_basis)
    end

    return RG
end

# Group ring with the basis the whole ball
function group_ring(G, half_radius::Integer, star_multiplication = false)
    S = Groups.gens(G)
    S = unique([S; inv.(S)])
    ball, _ = Groups.wlmetric_ball(S, one(G), radius = half_radius)
    return group_ring(G, ball; star_multiplication=star_multiplication, star_closed=true)
end
