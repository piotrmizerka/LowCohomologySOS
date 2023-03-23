using Revise
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
using JuMP
using SparseArrays
include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

sl3 = MatrixGroups.SpecialLinearGroup{3}(Int8)
F_sl_3_z = FreeGroup(alphabet(sl3))
e12, e13, e21, e23, e31, e32 = Groups.gens(F_sl_3_z)
quotient_hom = let source = F_sl_3_z, target = sl3
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end

relations = [
        e12 * e13 * e12^(-1) * e13^(-1),
        e13 * e12 * e13^(-1) * e12^(-1),
        e21 * e23 * e21^(-1) * e23^(-1),
        e23 * e21 * e23^(-1) * e21^(-1),
        e31 * e32 * e31^(-1) * e32^(-1),
        e32 * e31 * e32^(-1) * e31^(-1),
        e21 * e31 * e21^(-1) * e31^(-1),
        e31 * e21 * e31^(-1) * e21^(-1),
        e12 * e32 * e12^(-1) * e32^(-1),
        e32 * e12 * e32^(-1) * e12^(-1),
        e13 * e23 * e13^(-1) * e23^(-1),
        e23 * e13 * e23^(-1) * e13^(-1),
        e12 * e23 * e12^(-1) * e23^(-1) * e13^(-1),
        e13 * e32 * e13^(-1) * e32^(-1) * e12^(-1),
        e21 * e13 * e21^(-1) * e13^(-1) * e23^(-1),
        e23 * e31 * e23^(-1) * e31^(-1) * e21^(-1),
        e31 * e12 * e31^(-1) * e12^(-1) * e32^(-1),
        e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1),
]
d₁ = LowCohomologySOS.jacobian_matrix(relations, gens(F_sl_3_z))

S = let s = gens(sl3)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S, radius = 2)
RG = LowCohomologySOS.group_ring(sl3, half_basis, star_multiplication = false)
d0x = LowCohomologySOS.embed.(Ref(quotient_hom), LowCohomologySOS.d₀(parent(first(d₁)), gens(F_sl_3_z)), Ref(RG))
d1x = LowCohomologySOS.embed.(Ref(quotient_hom), d₁, Ref(RG))
a = d0x[1]
half_basis_restr = [one(sl3)]
for j in eachindex(d0x)
    for i in SparseArrays.nonzeroinds(d0x[j].coeffs)
        push!(half_basis_restr, RG.basis[i])
    end
end
for j in eachindex(d1x)
    @info j
    for i in SparseArrays.nonzeroinds(d1x[j].coeffs)
        push!(half_basis_restr, RG.basis[i])
    end
end
half_basis_restr = unique!([half_basis_restr; inv.(half_basis_restr)])
RGx = LowCohomologySOS.group_ring(sl3, half_basis_restr, star_multiplication = false)
d0xx = LowCohomologySOS.embed.(identity, d0x, Ref(RGx))
d1xx = LowCohomologySOS.embed.(identity, d1x, Ref(RGx))
Δ1⁺x = d1xx' * d1xx
Δ1⁻x = d0xx * d0xx'
RG_star = LowCohomologySOS.group_ring(sl3, half_basis_restr, star_multiplication = true)

Δ1⁺ = LowCohomologySOS.embed.(identity, Δ1⁺x, Ref(RG_star))
Δ1⁻ = LowCohomologySOS.embed.(identity, Δ1⁻x, Ref(RG_star))
sq, adj, op = LowCohomologySOS.sq_adj_op(Δ1⁻,gens(sl3))
Adj = 10*Δ1⁺+adj
I = [i ≠ j ? zero(parent(first(Adj))) : one(parent(first(Adj))) for i in 1:length(d0x), j in 1:length(d0x)]

sos_problem = LowCohomologySOS.sos_problem(Adj, I)
sos_problem = LowCohomologySOS.sos_problem(Δ1⁺+Δ1⁻, I)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 100_000))
JuMP.optimize!(sos_problem)
λ, Q = LowCohomologySOS.get_solution(sos_problem)

LowCohomologySOS.certify_sos_decomposition(Adj, I, λ, Q, half_basis_restr)
LowCohomologySOS.certify_sos_decomposition(Δ1⁺+Δ1⁻, I, λ, Q, half_basis_restr)

SL₃ℤ_data = (
    M = Δ1⁺+Δ1⁻,#Adj,
    order_unit = I,
    half_basis = half_basis_restr
)

solve_in_loop(
    sos_problem,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = SL₃ℤ_data
)
