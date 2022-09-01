using Groups

A = Alphabet([:x, :X, :y, :Y, :z, :Z], [2, 1, 4, 3, 6, 5])
F₃ = FreeGroup(A)
x, y, z = Groups.gens(F₃)

Π₀₁₋₁₋₁ = FPGroup(F₃, [x*y => y*x, x*z => z^(-1)*x, y*z => z^(-1)*y], maxrules=24)

S_Π₀₁₋₁₋₁ = let s = gens(Π₀₁₋₁₋₁)
    [s; inv.(s)]
end

const half_radius_Π₀₁₋₁₋₁ = 7
half_basis_Π₀₁₋₁₋₁, sizes_Π₀₁₋₁₋₁ = Groups.wlmetric_ball(S_Π₀₁₋₁₋₁, radius = half_radius_Π₀₁₋₁₋₁)

Π₁₁₁₁ = FPGroup(F₃, [x*y => z*y*x, x*z => z*x, y*z => z*y], maxrules=1000)

S_Π₁₁₁₁ = let s = gens(Π₁₁₁₁)
    [s; inv.(s)]
end

const half_radius_Π₁₁₁₁ = 7
half_basis_Π₁₁₁₁, sizes_Π₁₁₁₁ = Groups.wlmetric_ball(S_Π₁₁₁₁, radius = half_radius_Π₁₁₁₁)