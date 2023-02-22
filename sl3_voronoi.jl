using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "./")))
using Groups
using LowCohomologySOS
using JuMP
include(joinpath(@__DIR__, "./scripts/optimizers.jl"))
include(joinpath(@__DIR__, "./scripts/utils.jl"))

sl3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
S = gens(sl3)
S_inv =[S; inv.(S)]
basis, sizes = Groups.wlmetric_ball(S_inv, radius = 4)
sizes
half_basis = basis[1:sizes[4]]
RG = LowCohomologySOS.group_ring(sl3, half_basis, star_multiplication = false)

function gelt_from_matrix(M::AbstractMatrix, group_subset)
    for g in group_subset
        if MatrixGroups.matrix_repr(g) == M
            return g
        end
    end
end

function rg_elt(coeff_matrix_list)
    return sum(cm[1]*RG(gelt_from_matrix(cm[2],basis)) for cm in coeff_matrix_list)
end

d3₁₁ = rg_elt([
    (1,[1 1 0;0 1 0;0 0 1]),
    (-1,[1 0 0;1 1 0;0 0 1]),
    (-1,[1 0 0;0 1 0;0 0 1])
])
d3₁₂ = rg_elt([
    (1,[1 1 1; 0 1 0; 0 0 1]),
    (-1,[1 0 0; 1 1 1; 0 0 1]),
    (1,[1 0 0; 0 1 0; 1 1 1]),
    (-1,[1 0 0;0 1 0;0 0 1])
])
d4₁₁ = rg_elt([
    (-1,[1 0 0; 0 0 1; -1 -1 0]),
    (1,[1 0 0; 0 1 0; 1 0 1]),
    (-1,[1 0 0; 0 0 1; 0 -1 0]),
    (1,[1 0 0;0 1 0;0 0 1])
])
d4₂₁ = rg_elt([
    (1,[1 1 0; 0 0 1; 0 -1 0])
])

d3 = [d3₁₁ d3₁₂]
d4 = [
    d4₁₁;
    d4₂₁
]
d4*d4'
d3'*d3
Δ₃ = d4*d4'+d3'*d3
I₂ = [one(RG) zero(RG);zero(RG) one(RG)]

half_basis_prime = basis[1:sizes[2]]
RG_prime = LowCohomologySOS.group_ring(sl3, half_basis_prime, star_multiplication = true)
Δ₃ = LowCohomologySOS.embed.(identity, Δ₃, Ref(RG_prime))
I₂ = LowCohomologySOS.embed.(identity, I₂, Ref(RG_prime))

sl3_voronoi_data = (
    M = Δ₃,
    order_unit = I₂,
    half_basis = half_basis_prime
)

Δ₃_sos_problem = LowCohomologySOS.sos_problem(Δ₃, I₂)

solve_in_loop(
    Δ₃_sos_problem,
    logdir = "./logs_polaris_sl3_voronoi",
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = sl3_voronoi_data
)

