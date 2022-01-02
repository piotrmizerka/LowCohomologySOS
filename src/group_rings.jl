Base.:/(G::Groups.AbstractFPGroup, rels::AbstractVector{<:FPGroupElement}) = FPGroup(G, [r=>one(G) for r in rels])
Base.adjoint(X::AlgebraElement) = StarAlgebras.star(X)
Base.copy(X::AlgebraElement) = AlgebraElement(copy(StarAlgebras.coeffs(X)),parent(X))
StarAlgebras.star(g::Groups.GroupElement) = inv(g)

# Group ring with the basis given by the prescribed support
function group_ring(G, half_basis, star_multiplication = false)
    star_closed_support = unique!([half_basis; star.(half_basis)])
    basis_elements = let basis = star_closed_support, f = star_multiplication ? star : identity
        unique!([basis;[f(a)*b for a in basis for b in basis]])
    end
    group_ring_basis = StarAlgebras.Basis{Int}(basis_elements)
    tmstr = StarAlgebras.CachedMTable{star_multiplication}(group_ring_basis, table_size = (length(star_closed_support), length(star_closed_support)))
    RG = StarAlgebra(G, group_ring_basis, tmstr)

    return RG
end

# Group ring with the basis the whole ball
function group_ring(G, half_radius::Integer, star_multiplication = false)
    S = Groups.gens(G)
    S = unique([S; inv.(S)])
    ball, sizes = Groups.wlmetric_ball(S, one(G), radius = 2*half_radius)
    b = StarAlgebras.Basis{Int}(ball)
    tmstr = StarAlgebras.CachedMTable{star_multiplication}(b, table_size = (sizes[half_radius], sizes[half_radius]))
    RG = StarAlgebra(G, b, tmstr)

    return RG
end
