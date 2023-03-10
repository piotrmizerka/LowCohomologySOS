include(joinpath(@__DIR__, "voronoi_utils.jl"));

S_m31, S_m32, S_m4 = saturate(S,m31_stab), saturate(S,m32_stab), saturate(S,m4_stab);
S_sat = collect(union(S_m31, S_m32, S_m4,m31_stab, m32_stab, m4_stab));
S_sat_inv = collect(Set([S_sat; inv.(S_sat)]));
ball2_sat, sizes_sat = Groups.wlmetric_ball(S_sat_inv, radius = 2);
outsiders_d4₁₁ = [gelt_from_matrix(M,ball8) for M in [
    [-1 0 0; 1 0 1; 1 1 0],[0 0 1; -1 0 -1; -1 -1 0],
    [-1 0 -1; 1 0 0; -1 -1 0],[1 0 1; 0 0 -1; 1 1 0]
]]
# half_basis_sat = collect(union(ball2_sat, outsiders_d4₁₁))
temp_basis = collect(union(ball2_sat, outsiders_d4₁₁))
d411_array = [
    [1 0 0; 0 0 1; -1 -1 0], [-1 0 0; 0 0 -1; -1 -1 0], [-1 0 0; 1 0 1; 1 1 0], 
    [1 0 0; -1 0 -1; 1 1 0], [0 0 1; 1 0 0; 1 1 0], [0 0 -1; -1 0 0; 1 1 0], 
    [0 0 -1; 1 0 1; -1 -1 0], [0 0 1; -1 0 -1; -1 -1 0], [-1 0 -1; 1 0 0; -1 -1 0], 
    [1 0 1; -1 0 0; -1 -1 0], [-1 0 -1; 0 0 1; 1 1 0], [1 0 1; 0 0 -1; 1 1 0]
]
d412_array = [
    [1 0 0; 0 1 0; 1 0 1], [-1 0 0; 0 -1 0; 1 0 1], [-1 0 0; 1 1 0; -1 0 -1], 
    [1 0 0; -1 -1 0; -1 0 -1], [0 1 0; 1 0 0; -1 0 -1], [0 -1 0; -1 0 0; -1 0 -1], 
    [0 -1 0; 1 1 0; 1 0 1], [0 1 0; -1 -1 0; 1 0 1], [-1 -1 0; 1 0 0; 1 0 1], 
    [1 1 0; -1 0 0; 1 0 1], [-1 -1 0; 0 1 0; -1 0 -1], [1 1 0; 0 -1 0; -1 0 -1]
]
d413_array = [
    [1 0 0; 0 0 1; 0 -1 0], [-1 0 0; 0 0 -1; 0 -1 0], [-1 0 0; 1 0 1; 0 1 0], 
    [1 0 0; -1 0 -1; 0 1 0], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; -1 0 0; 0 1 0], 
    [0 0 -1; 1 0 1; 0 -1 0], [0 0 1; -1 0 -1; 0 -1 0], [-1 0 -1; 1 0 0; 0 -1 0], 
    [1 0 1; -1 0 0; 0 -1 0], [-1 0 -1; 0 0 1; 0 1 0], [1 0 1; 0 0 -1; 0 1 0]
]
d414_array = [
    [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 1 1 0; 0 0 -1], 
    [1 0 0; -1 -1 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 -1], 
    [0 -1 0; 1 1 0; 0 0 1], [0 1 0; -1 -1 0; 0 0 1], [-1 -1 0; 1 0 0; 0 0 1], 
    [1 1 0; -1 0 0; 0 0 1], [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; 0 0 -1]
]
d42_array = [
    [0 1 0; 0 0 -1; -1 -1 0], [0 -1 0; 0 0 1; -1 0 -1], [0 -1 0; 1 1 0; 0 0 1], 
    [0 1 0; -1 -1 0; 1 0 1], [0 1 0; 1 0 1; 0 0 -1], [0 -1 0; -1 0 -1; 1 1 0], 
    [0 0 1; 0 -1 0; 1 1 0], [0 0 -1; 0 1 0; 1 0 1], [0 0 -1; -1 -1 0; 0 1 0], 
    [0 0 1; 1 1 0; -1 0 -1], [0 0 1; -1 0 -1; 0 -1 0], [0 0 -1; 1 0 1; -1 -1 0], 
    [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; -1 0 -1], [1 1 0; 0 0 1; 0 -1 0], 
    [-1 -1 0; 0 0 -1; 1 0 1], [-1 -1 0; 1 0 1; 0 1 0], [1 1 0; -1 0 -1; 0 0 1], 
    [-1 0 -1; 0 -1 0; 0 0 1], [1 0 1; 0 1 0; -1 -1 0], [1 0 1; 0 0 -1; 0 1 0], 
    [-1 0 -1; 0 0 1; 1 1 0], [-1 0 -1; 1 1 0; 0 -1 0], [1 0 1; -1 -1 0; 0 0 -1]
]
d51_array = [
    [1 1 1; 0 -1 0; 0 0 -1], [-1 -1 -1; 0 1 0; 1 1 0], [-1 -1 -1; 0 0 1; 0 1 0], [1 1 1; 0 0 -1; -1 0 -1], 
    [1 1 1; -1 -1 0; 0 -1 0], [-1 -1 -1; 1 1 0; 1 0 1], [-1 -1 -1; 1 0 1; 0 0 1], [1 1 1; -1 0 -1; -1 -1 0]
]
d52_array = [
    [0 0 1; 1 0 0; 1 1 0], [0 0 -1; -1 0 0; 1 1 1], [0 0 -1; -1 -1 0; -1 0 0], [0 0 1; 1 1 0; -1 0 -1], 
    [0 0 -1; 1 0 1; -1 -1 0], [0 0 1; -1 0 -1; -1 -1 -1], [0 0 1; -1 -1 -1; 1 0 0], [0 0 -1; 1 1 1; 1 0 1]
]
d53_array = [
    [0 -1 0; -1 0 0; -1 0 -1], [0 1 0; 1 0 0; -1 -1 -1], [0 1 0; -1 -1 0; 1 0 1], [0 -1 0; 1 1 0; 1 1 1], 
    [0 1 0; 1 0 1; 1 0 0], [0 -1 0; -1 0 -1; 1 1 0], [0 -1 0; 1 1 1; -1 0 0], [0 1 0; -1 -1 -1; -1 -1 0]
]
d54_array = [
    [-1 0 -1; 1 0 0; 0 -1 0], [1 0 1; -1 0 0; -1 -1 -1], [1 0 1; 0 1 0; -1 0 0], [-1 0 -1; 0 -1 0; 0 0 1], 
    [1 0 1; 0 0 -1; 0 1 0], [-1 0 -1; 0 0 1; 1 1 1], [-1 0 -1; 1 1 1; 1 0 0], [1 0 1; -1 -1 -1; 0 0 -1]
]
d55_array = [
    [1 1 0; -1 0 0; 0 0 1], [-1 -1 0; 1 0 0; 1 1 1], [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; -1 -1 -1], 
    [-1 -1 0; 0 0 -1; 1 0 0], [1 1 0; 0 0 1; 0 -1 0], [1 1 0; -1 -1 -1; -1 0 0], [-1 -1 0; 1 1 1; 0 1 0]
]
d56_array = [
    [1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 1 0 1], [-1 0 0; 0 0 -1; 0 -1 0], [1 0 0; 0 0 1; -1 -1 0], 
    [-1 0 0; 1 1 0; 0 0 -1], [1 0 0; -1 -1 0; -1 0 -1], [1 0 0; -1 0 -1; 0 1 0], [-1 0 0; 1 0 1; 1 1 0]
]
d411 = [gelt_from_matrix(M,temp_basis) for M in d411_array]
d412 = [gelt_from_matrix(M,temp_basis) for M in d412_array]
d413 = [gelt_from_matrix(M,temp_basis) for M in d413_array]
d414 = [gelt_from_matrix(M,temp_basis) for M in d414_array]
d42 = [gelt_from_matrix(M,temp_basis) for M in d42_array]
d51 = [gelt_from_matrix(M,temp_basis) for M in d51_array]
d52 = [gelt_from_matrix(M,temp_basis) for M in d52_array]
d53 = [gelt_from_matrix(M,temp_basis) for M in d53_array]
d54 = [gelt_from_matrix(M,temp_basis) for M in d54_array]
d55 = [gelt_from_matrix(M,temp_basis) for M in d55_array]
d56 = [gelt_from_matrix(M,temp_basis) for M in d56_array]

half_basis_sat = unique(union(d411, d412, d413, d414, d42, d51, d52, d53, d54, d55, d56, m4_stab))
half_basis_sat = unique([half_basis_sat;inv.(half_basis_sat)])
RG = LowCohomologySOS.group_ring(sl3, half_basis_sat, star_multiplication = false)

d4₁ = rg_elt([
    (-1,averaged_rep(d411_array, half_basis_sat, RG)),
    (1,averaged_rep(d412_array, half_basis_sat, RG)),
    (-1,averaged_rep(d413_array, half_basis_sat, RG)),
    (1,averaged_rep(d414_array, half_basis_sat, RG))
])
d4₂ = rg_elt([
    (1,averaged_rep(d42_array, half_basis_sat, RG))
])
d5 = rg_elt([
    (1,averaged_rep(d51_array, half_basis_sat, RG)),
    (-1,averaged_rep(d52_array, half_basis_sat, RG)),
    (1,averaged_rep(d53_array, half_basis_sat, RG)),
    (-1,averaged_rep(d54_array, half_basis_sat, RG)),
    (1,averaged_rep(d55_array, half_basis_sat, RG)),
    (-1,averaged_rep(d56_array, half_basis_sat, RG))
])

d₄ = [
    d4₁;
    d4₂
]
d₅ = reshape([d5], 1, 1)
m4_stab_part = reshape([one(RG)-averaged_rep(m4_arrays, half_basis_sat, RG)], 1, 1)
Δ₄x = reshape([d₄'*d₄],1,1)+d₅*d₅'+m4_stab_part
RG_star = LowCohomologySOS.group_ring(sl3, half_basis_sat, star_multiplication = true)
Δ₄ = LowCohomologySOS.embed.(identity, Δ₄x, Ref(RG_star))
I = reshape([one(RG_star)], 1, 1)

sos_problem = LowCohomologySOS.sos_problem(Δ₄, I, 1.93)

JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 4_000))
JuMP.optimize!(sos_problem)

sl3_voronoi_data = (
    M = Δ₄,
    order_unit = I,
    half_basis = half_basis_sat
)
solve_in_loop(
    sos_problem,
    logdir = joinpath(@__DIR__, "logs_voronoi_delta_4"),
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = sl3_voronoi_data
)