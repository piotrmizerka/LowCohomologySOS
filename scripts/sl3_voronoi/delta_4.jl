include(joinpath(@__DIR__, "voronoi_utils.jl"));

S_m31, S_m32, S_m4 = saturate(S,m31_stab), saturate(S,m32_stab), saturate(S,m4_stab);
S_sat = collect(union(S_m31, S_m32, S_m4,m31_stab, m32_stab, m4_stab));
S_sat_inv = collect(Set([S_sat; inv.(S_sat)]));
ball2_sat, sizes_sat = Groups.wlmetric_ball(S_sat_inv, radius = 2);
outsiders_d4₁₁ = [gelt_from_matrix(M,ball8) for M in [
    [-1 0 0; 1 0 1; 1 1 0],[0 0 1; -1 0 -1; -1 -1 0],
    [-1 0 -1; 1 0 0; -1 -1 0],[1 0 1; 0 0 -1; 1 1 0]
]]
half_basis_sat = collect(union(ball2_sat, outsiders_d4₁₁))
RG = LowCohomologySOS.group_ring(sl3, half_basis_sat, star_multiplication = false)

d4₁ = rg_elt([
    (-1,averaged_rep(
        [[1 0 0; 0 0 1; -1 -1 0], [-1 0 0; 0 0 -1; -1 -1 0], [-1 0 0; 1 0 1; 1 1 0], 
        [1 0 0; -1 0 -1; 1 1 0], [0 0 1; 1 0 0; 1 1 0], [0 0 -1; -1 0 0; 1 1 0], 
        [0 0 -1; 1 0 1; -1 -1 0], [0 0 1; -1 0 -1; -1 -1 0], [-1 0 -1; 1 0 0; -1 -1 0], 
        [1 0 1; -1 0 0; -1 -1 0], [-1 0 -1; 0 0 1; 1 1 0], [1 0 1; 0 0 -1; 1 1 0]],
        half_basis_sat, RG
    )),
    (1,averaged_rep(
        [[1 0 0; 0 1 0; 1 0 1], [-1 0 0; 0 -1 0; 1 0 1], [-1 0 0; 1 1 0; -1 0 -1], 
        [1 0 0; -1 -1 0; -1 0 -1], [0 1 0; 1 0 0; -1 0 -1], [0 -1 0; -1 0 0; -1 0 -1], 
        [0 -1 0; 1 1 0; 1 0 1], [0 1 0; -1 -1 0; 1 0 1], [-1 -1 0; 1 0 0; 1 0 1], 
        [1 1 0; -1 0 0; 1 0 1], [-1 -1 0; 0 1 0; -1 0 -1], [1 1 0; 0 -1 0; -1 0 -1]],
        half_basis_sat, RG

    )),
    (-1,averaged_rep(
        [[1 0 0; 0 0 1; 0 -1 0], [-1 0 0; 0 0 -1; 0 -1 0], [-1 0 0; 1 0 1; 0 1 0], 
        [1 0 0; -1 0 -1; 0 1 0], [0 0 1; 1 0 0; 0 1 0], [0 0 -1; -1 0 0; 0 1 0], 
        [0 0 -1; 1 0 1; 0 -1 0], [0 0 1; -1 0 -1; 0 -1 0], [-1 0 -1; 1 0 0; 0 -1 0], 
        [1 0 1; -1 0 0; 0 -1 0], [-1 0 -1; 0 0 1; 0 1 0], [1 0 1; 0 0 -1; 0 1 0]],
        half_basis_sat, RG
    )),
    (1,averaged_rep(
        [[1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 0 0 1], [-1 0 0; 1 1 0; 0 0 -1], 
        [1 0 0; -1 -1 0; 0 0 -1], [0 1 0; 1 0 0; 0 0 -1], [0 -1 0; -1 0 0; 0 0 -1], 
        [0 -1 0; 1 1 0; 0 0 1], [0 1 0; -1 -1 0; 0 0 1], [-1 -1 0; 1 0 0; 0 0 1], 
        [1 1 0; -1 0 0; 0 0 1], [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; 0 0 -1]],
        half_basis_sat, RG
    ))
])
d4₂ = rg_elt([
    (1,averaged_rep(
        [[0 1 0; 0 0 -1; -1 -1 0], [0 -1 0; 0 0 1; -1 0 -1], [0 -1 0; 1 1 0; 0 0 1], 
        [0 1 0; -1 -1 0; 1 0 1], [0 1 0; 1 0 1; 0 0 -1], [0 -1 0; -1 0 -1; 1 1 0], 
        [0 0 1; 0 -1 0; 1 1 0], [0 0 -1; 0 1 0; 1 0 1], [0 0 -1; -1 -1 0; 0 1 0], 
        [0 0 1; 1 1 0; -1 0 -1], [0 0 1; -1 0 -1; 0 -1 0], [0 0 -1; 1 0 1; -1 -1 0], 
        [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; -1 0 -1], [1 1 0; 0 0 1; 0 -1 0], 
        [-1 -1 0; 0 0 -1; 1 0 1], [-1 -1 0; 1 0 1; 0 1 0], [1 1 0; -1 0 -1; 0 0 1], 
        [-1 0 -1; 0 -1 0; 0 0 1], [1 0 1; 0 1 0; -1 -1 0], [1 0 1; 0 0 -1; 0 1 0], 
        [-1 0 -1; 0 0 1; 1 1 0], [-1 0 -1; 1 1 0; 0 -1 0], [1 0 1; -1 -1 0; 0 0 -1]],
        half_basis_sat, RG
    ))
])
d5 = rg_elt([
    (1,averaged_rep(
        [[1 1 1; 0 -1 0; 0 0 -1], [-1 -1 -1; 0 1 0; 1 1 0], [-1 -1 -1; 0 0 1; 0 1 0], [1 1 1; 0 0 -1; -1 0 -1], 
        [1 1 1; -1 -1 0; 0 -1 0], [-1 -1 -1; 1 1 0; 1 0 1], [-1 -1 -1; 1 0 1; 0 0 1], [1 1 1; -1 0 -1; -1 -1 0]],
        half_basis_sat, RG
    )),
    (-1,averaged_rep(
        [[0 0 1; 1 0 0; 1 1 0], [0 0 -1; -1 0 0; 1 1 1], [0 0 -1; -1 -1 0; -1 0 0], [0 0 1; 1 1 0; -1 0 -1], 
        [0 0 -1; 1 0 1; -1 -1 0], [0 0 1; -1 0 -1; -1 -1 -1], [0 0 1; -1 -1 -1; 1 0 0], [0 0 -1; 1 1 1; 1 0 1]],
        half_basis_sat, RG
    )),
    (1,averaged_rep(
        [[0 -1 0; -1 0 0; -1 0 -1], [0 1 0; 1 0 0; -1 -1 -1], [0 1 0; -1 -1 0; 1 0 1], [0 -1 0; 1 1 0; 1 1 1], 
        [0 1 0; 1 0 1; 1 0 0], [0 -1 0; -1 0 -1; 1 1 0], [0 -1 0; 1 1 1; -1 0 0], [0 1 0; -1 -1 -1; -1 -1 0]],
        half_basis_sat, RG
    )),
    (-1,averaged_rep(
        [[-1 0 -1; 1 0 0; 0 -1 0], [1 0 1; -1 0 0; -1 -1 -1], [1 0 1; 0 1 0; -1 0 0], [-1 0 -1; 0 -1 0; 0 0 1], 
        [1 0 1; 0 0 -1; 0 1 0], [-1 0 -1; 0 0 1; 1 1 1], [-1 0 -1; 1 1 1; 1 0 0], [1 0 1; -1 -1 -1; 0 0 -1]],
        half_basis_sat, RG
    )),
    (1,averaged_rep(
        [[1 1 0; -1 0 0; 0 0 1], [-1 -1 0; 1 0 0; 1 1 1], [-1 -1 0; 0 1 0; 0 0 -1], [1 1 0; 0 -1 0; -1 -1 -1], 
        [-1 -1 0; 0 0 -1; 1 0 0], [1 1 0; 0 0 1; 0 -1 0], [1 1 0; -1 -1 -1; -1 0 0], [-1 -1 0; 1 1 1; 0 1 0]],
        half_basis_sat, RG
    )),
    (-1,averaged_rep(
        [[1 0 0; 0 1 0; 0 0 1], [-1 0 0; 0 -1 0; 1 0 1], [-1 0 0; 0 0 -1; 0 -1 0], [1 0 0; 0 0 1; -1 -1 0], 
        [-1 0 0; 1 1 0; 0 0 -1], [1 0 0; -1 -1 0; -1 0 -1], [1 0 0; -1 0 -1; 0 1 0], [-1 0 0; 1 0 1; 1 1 0]],
        half_basis_sat, RG
    ))
])

m31_stab_part = one(RG)-averaged_rep(m31_arrays, half_basis_sat, RG)
m32_stab_part = one(RG)-averaged_rep(m32_arrays, half_basis_sat, RG)
m4_stab_part = one(RG)-averaged_rep(m4_arrays, half_basis_sat, RG)
d₄ = [d4₁+m31_stab_part d4₂+m32_stab_part]
d₅ = reshape([d5+m4_stab_part], 1, 1)

Δ₄x = d₅'*d₅+d₄*d₄'
RG_star = LowCohomologySOS.group_ring(sl3, half_basis_sat, star_multiplication = true)
Δ₄ = LowCohomologySOS.embed.(identity, Δ₄x, Ref(RG_star))
I = reshape([one(RG_star)], 1, 1)

sos_problem = LowCohomologySOS.sos_problem(Δ₄, I)

sl3_voronoi_data = (
    M = Δ₄,
    order_unit = I,
    half_basis = half_basis_sat
)
solve_in_loop(
    sos_problem,
    logdir = joinpath(@__DIR__, "logs_voronoi"),
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = sl3_voronoi_data
)