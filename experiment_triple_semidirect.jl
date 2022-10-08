using Revise
using Groups
using LowCohomologySOS

include(joinpath(@__DIR__, "scripts", "optimizers.jl"))
include(joinpath(@__DIR__, "scripts", "utils.jl"))

function check_vanishing(G::Groups.AbstractFPGroup, d₁, d₂, 
    half_basis)
    ℝG_twisted = LowCohomologySOS.group_ring(G, half_basis, star_multiplication = true)

    Δ₂⁺ = d₂'*d₂
    Δ₂⁻ = d₁*d₁'

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

# Π₀₋₁₁₋₁ ##############################################################################################
G = 𝔹₃ = Π₀₋₁₁₋₁ = FPGroup(F₃, [x*y => y^(-1)*x, x*z => z*x, y*z => z^(-1)*y])

S = let s = gens(G)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)
ℝG = LowCohomologySOS.group_ring(G, half_basis)
quotient_hom_G = quotient_hom(F₃, G)
xx, yy, zz = quotient_hom_G.([x, y, z])

d₁ = [
        ℝG(yy^(-1))-one(ℝG) ℝG((xx*yy)^(-1))+one(ℝG) zero(ℝG);
        ℝG(zz^(-1))-one(ℝG) zero(ℝG) one(ℝG)-ℝG(xx^(-1));
        zero(ℝG) ℝG(zz^(-1))-one(ℝG) ℝG((yy*zz)^(-1))+one(ℝG)
]

d₂ = [
        ℝG(zz^(-1))-one(ℝG) ℝG((yy*zz)^(-1))+one(ℝG) ℝG((xx*yy*zz)^(-1))-one(ℝG)
]

check_vanishing(G, d₁, d₂, half_basis)


# Π₁₋₁₁₋₁ ##############################################################################################
G = 𝔹₄ = Π₁₋₁₁₋₁ = FPGroup(F₃, [x*y => z*y^(-1)*x, x*z => z*x, y*z => z^(-1)*y])

S = let s = gens(G)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)
ℝG = LowCohomologySOS.group_ring(G, half_basis)
quotient_hom_G = quotient_hom(F₃, G)
xx, yy, zz = quotient_hom_G.([x, y, z])

d₁ = [
        ℝG((yy*zz)^(-1))-one(ℝG) ℝG((xx*yy*zz)^(-1))+ℝG(zz^(-1)) one(ℝG);
        ℝG(zz^(-1))-one(ℝG) zero(ℝG) one(ℝG)-ℝG(xx^(-1));
        zero(ℝG) ℝG(zz^(-1))-one(ℝG) ℝG((yy*zz)^(-1))+one(ℝG)
]

d₂ = [
        ℝG(zz^(-1))-one(ℝG) ℝG((yy*zz^2)^(-1))+one(ℝG) ℝG((xx*yy*zz^2)^(-1))-ℝG(zz^(-1))
]

check_vanishing(G, d₁, d₂, half_basis)
