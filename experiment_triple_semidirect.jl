using Revise
using Groups
using LowCohomologySOS

include(joinpath(@__DIR__, "scripts", "optimizers.jl"))
include(joinpath(@__DIR__, "scripts", "utils.jl"))

function check_vanishing(G::Groups.AbstractFPGroup, d₁, d₂, 
    half_basis)
    ℝG_twisted = LowCohomologySOS.group_ring(G, half_basis, star_multiplication = true)

    # right multiplication convention reverses order of stars in multiplication ??
    Δ₂⁺ = d₂*d₂'
    Δ₂⁻ = d₁'*d₁

    Δ₂ = Δ₂⁺+Δ₂⁻
    Δ₂_twisted = LowCohomologySOS.embed.(identity, Δ₂, Ref(ℝG_twisted))
    I₃ = [i ≠ j ? zero(ℝG_twisted) : one(ℝG_twisted) for i in 1:3, j in 1:3]

    Δ₂_sos_problem = LowCohomologySOS.sos_problem_matrix(
        Δ₂_twisted,
        I₃
    )
    G_data = let
        (
            M = Δ₂_twisted,
            order_unit = I₃,
            half_basis = half_basis,
            RG = parent(first(Δ₂_twisted)),
        )
    end

    solve_in_loop(
        Δ₂_sos_problem,
        logdir = "./logs",
        optimizer = scs_opt(eps = 1e-5, max_iters = 20_000),
        data = G_data
    )
end

function quotient_hom(source, target)
    result = let source = source, target = target
        function f(i, source, target)
            if alphabet(source) == alphabet(target)
                Groups.word_type(target)([i])
            else
                throw("Unsupported")
            end
        end
        Groups.Homomorphism(f, source, target)
    end
    return result
end

A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
F₃ = FreeGroup(A)
x, y, z = Groups.gens(F₃)

const half_radius = 5


# Π₀₁₋₁₋₁ ##############################################################################################
G = Π₀₁₋₁₋₁ = FPGroup(F₃, [x*y => y*x, x*z => z^(-1)*x, y*z => z^(-1)*y], maxrules=24)

S = let s = gens(G)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)
ℝG = LowCohomologySOS.group_ring(G, half_basis)
quotient_hom_G = quotient_hom(F₃, G)
xx, yy, zz = quotient_hom_G.([x, y, z])
# Is dᵢ=∂ᵢ*???
d₁ = [
        one(ℝG)-ℝG(yy^(-1)) ℝG(zz^(-1))-one(ℝG) zero(ℝG);
        ℝG(xx^(-1))-one(ℝG) zero(ℝG) ℝG(zz^(-1))-one(ℝG);
        zero(ℝG) one(ℝG)+ℝG(xx^(-1)*zz^(-1)) one(ℝG)+ℝG(yy^(-1)*zz^(-1))
]

d₂ = [
        one(ℝG)-ℝG(zz^(-1));
        one(ℝG)+ℝG(yy^(-1)*zz^(-1));
        -ℝG(xx^(-1)*zz^(-1))-one(ℝG)
]

check_vanishing(G, d₁, d₂, half_basis)


# Π₀₋₁₁₋₁ ##############################################################################################
G = Π₀₋₁₁₋₁ = FPGroup(F₃, [x*y => y^(-1)*x, x*z => z*x, y*z => z^(-1)*y])

S = let s = gens(G)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)
ℝG = LowCohomologySOS.group_ring(G, half_basis)
quotient_hom_G = quotient_hom(F₃, G)
xx, yy, zz = quotient_hom_G.([x, y, z])

d₁ = [
        -one(ℝG)+ℝG(yy^(-1)) -ℝG(zz^(-1))+one(ℝG) zero(ℝG);
        one(ℝG)+ℝG(xx^(-1)*yy^(-1)) zero(ℝG) ℝG(zz^(-1))-one(ℝG);
        zero(ℝG) -one(ℝG)+ℝG(xx^(-1)) one(ℝG)+ℝG(yy^(-1)*zz^(-1))
]

d₂ = [
        one(ℝG)-ℝG(zz^(-1));
        one(ℝG)+ℝG(yy^(-1)*zz^(-1));
        ℝG(zz^(-1))-ℝG(xx^(-1)*yy^(-1)*zz^(-1))
]

check_vanishing(G, d₁, d₂, half_basis)


# Π₁₋₁₁₋₁ ##############################################################################################
G = Π₁₋₁₁₋₁ = FPGroup(F₃, [x*y => z*y^(-1)*x, x*z => z*x, y*z => z^(-1)*y])

S = let s = gens(G)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)
ℝG = LowCohomologySOS.group_ring(G, half_basis)
quotient_hom_G = quotient_hom(F₃, G)
xx, yy, zz = quotient_hom_G.([x, y, z])

d₁ = [
        -one(ℝG)+ℝG(yy^(-1)*zz^(-1)) one(ℝG)-ℝG(zz^(-1)) zero(ℝG);
        ℝG(zz^(-1))+ℝG(xx^(-1)*yy^(-1)*zz^(-1)) zero(ℝG) ℝG(zz^(-1))-one(ℝG);
        one(ℝG) -one(ℝG)+ℝG(xx^(-1)) one(ℝG)+ℝG(yy^(-1)*zz^(-1))
]

d₂ = [
        one(ℝG)-ℝG(zz^(-1));
        one(ℝG)+ℝG(yy^(-1)*(zz^(-1))^2);
        ℝG(zz^(-1))-ℝG(xx^(-1)*yy^(-1)*(zz^(-1))^2)
]

check_vanishing(G, d₁, d₂, half_basis)


# check reducibility using lower Laplacian - caveat: the certification argument doees not work here!!! It has to be figured out!!!
Δ₂⁻_square = LowCohomologySOS.embed.(identity, Δ₂⁻^2, Ref(ℝΠ₀₁₋₁₋₁x))
Δ₂⁻x = LowCohomologySOS.embed.(identity, Δ₂⁻, Ref(ℝΠ₀₁₋₁₋₁x))

Π₀₁₋₁₋₁x_reducibility_lower_sos_problem = LowCohomologySOS.sos_problem_matrix(
        Δ₂⁻_square,
        Δ₂⁻x
)
Π₀₁₋₁₋₁_reducibility_lower_data = let
    (
        M = Δ₂⁻_square,
        order_unit = Δ₂⁻x,
        half_basis = half_basis_Π₀₁₋₁₋₁,
        RG = parent(first(Δ₂⁻x)),
    )
end
solve_in_loop(
    Π₀₁₋₁₋₁x_reducibility_lower_sos_problem,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-5, max_iters = 100_000),
    data = Π₀₁₋₁₋₁_reducibility_lower_data
)

# check reducibility using upper Laplacian TODO
Δ₁⁺ = d₁*d₁'
Δ₁⁺_square = LowCohomologySOS.embed.(identity, Δ₁⁺^2, Ref(ℝΠ₀₁₋₁₋₁x))
Δ₁⁺x = LowCohomologySOS.embed.(identity, Δ₁⁺, Ref(ℝΠ₀₁₋₁₋₁x))
