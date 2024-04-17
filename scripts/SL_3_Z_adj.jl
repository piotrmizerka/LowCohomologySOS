using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using JuMP
using LowCohomologySOS

include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

Adj, Iₙ, half_basis = let half_radius = 2
    SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    S = let s = gens(SL₃ℤ)
        [s; inv.(s)]
    end
    half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)

    F_sl_3_z = FreeGroup(alphabet(SL₃ℤ))
    e12, e13, e21, e23, e31, e32 = Groups.gens(F_sl_3_z)

    quotient_hom = let source = F_sl_3_z, target = SL₃ℤ
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

    Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)
    sq, adj, op = LowCohomologySOS.sq_adj_op(Δ₁⁻,gens(SL₃ℤ))
    Δ₁⁺+adj, Iₙ, half_basis
end

SL₃ℤ_data = (
    M = Adj,
    order_unit = Iₙ,
    half_basis = half_basis
)

sos_problem = LowCohomologySOS.sos_problem(Adj, Iₙ)

solve_in_loop(
    sos_problem,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 100_000),
    data = SL₃ℤ_data
)
