using Groups

A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
F₃ = FreeGroup(A)
x, y, z = Groups.gens(F₃)

# Π₀₁₋₁₋₁ ##############################################################################################
Π₀₁₋₁₋₁ = FPGroup(F₃, [x*y => y*x, x*z => z^(-1)*x, y*z => z^(-1)*y], maxrules=24)

S_Π₀₁₋₁₋₁ = let s = gens(Π₀₁₋₁₋₁)
    [s; inv.(s)]
end

const half_radius_Π₀₁₋₁₋₁ = 4
half_basis_Π₀₁₋₁₋₁, sizes_Π₀₁₋₁₋₁ = Groups.wlmetric_ball(S_Π₀₁₋₁₋₁, radius = half_radius_Π₀₁₋₁₋₁)

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

quotient_hom_Π₀₁₋₁₋₁ = quotient_hom(F₃, Π₀₁₋₁₋₁)
xx, yy, zz = quotient_hom_Π₀₁₋₁₋₁.([x, y, z])

using LowCohomologySOS

ℝΠ₀₁₋₁₋₁ = LowCohomologySOS.group_ring(Π₀₁₋₁₋₁, half_basis_Π₀₁₋₁₋₁)
x∈Π₀₁₋₁₋₁

# Is dᵢ=∂ᵢ*???
d₁ = [
        one(ℝΠ₀₁₋₁₋₁)-ℝΠ₀₁₋₁₋₁(yy^(-1)) ℝΠ₀₁₋₁₋₁(zz^(-1))-one(ℝΠ₀₁₋₁₋₁) zero(ℝΠ₀₁₋₁₋₁);
        ℝΠ₀₁₋₁₋₁(xx^(-1))-one(ℝΠ₀₁₋₁₋₁) zero(ℝΠ₀₁₋₁₋₁) ℝΠ₀₁₋₁₋₁(zz^(-1))-one(ℝΠ₀₁₋₁₋₁);
        zero(ℝΠ₀₁₋₁₋₁) one(ℝΠ₀₁₋₁₋₁)+ℝΠ₀₁₋₁₋₁(xx^(-1)*zz^(-1)) one(ℝΠ₀₁₋₁₋₁)+ℝΠ₀₁₋₁₋₁(yy^(-1)*zz^(-1))
]

d₂ = [
        one(ℝΠ₀₁₋₁₋₁)-ℝΠ₀₁₋₁₋₁(zz^(-1));
        one(ℝΠ₀₁₋₁₋₁)+ℝΠ₀₁₋₁₋₁(yy^(-1)*zz^(-1));
        -ℝΠ₀₁₋₁₋₁(xx^(-1)*zz^(-1))-one(ℝΠ₀₁₋₁₋₁)
]

# right multiplication convention reverses order of stars in multiplication ??
Δ₂⁺ = d₂*d₂'
Δ₂⁻ = d₁'*d₁

# check vanishing
Δ₂ = Δ₂⁺+Δ₂⁻

ℝΠ₀₁₋₁₋₁x = LowCohomologySOS.group_ring(Π₀₁₋₁₋₁, half_basis_Π₀₁₋₁₋₁, star_multiplication = true)
Δ₂x = LowCohomologySOS.embed.(identity, Δ₂, Ref(ℝΠ₀₁₋₁₋₁x))
I₃ = [i ≠ j ? zero(ℝΠ₀₁₋₁₋₁x) : one(ℝΠ₀₁₋₁₋₁x) for i in 1:3, j in 1:3]

Δ₂_sos_problem = LowCohomologySOS.sos_problem_matrix(
        Δ₂x,
        I₃
)
Π₀₁₋₁₋₁_data = let
    (
        M = Δ₂x,
        order_unit = I₃,
        half_basis = half_basis_Π₀₁₋₁₋₁,
        RG = parent(first(Δ₂x)),
    )
end

include(joinpath(@__DIR__, "scripts", "optimizers.jl"))
include(joinpath(@__DIR__, "scripts", "utils.jl"))

solve_in_loop(
    Δ₂_sos_problem,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-5, max_iters = 100_000),
    data = Π₀₁₋₁₋₁_data
)

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


# Π₁₁₁₁ ###### TODO - compactify the code first!
Π₁₁₁₁ = FPGroup(F₃, [x*y => z*y*x, x*z => z*x, y*z => z*y], maxrules=1000)

S_Π₁₁₁₁ = let s = gens(Π₁₁₁₁)
    [s; inv.(s)]
end

const half_radius_Π₁₁₁₁ = 7
half_basis_Π₁₁₁₁, sizes_Π₁₁₁₁ = Groups.wlmetric_ball(S_Π₁₁₁₁, radius = half_radius_Π₁₁₁₁)