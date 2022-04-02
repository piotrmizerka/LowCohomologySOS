using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
using PropertyT_new
using JuMP
include(joinpath(@__DIR__, "optimizers.jl"))

Δ₁x, Iₙ, half_basis = let half_radius = 2
    SL(n, R) = PropertyT_new.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    S = let s = gens(SL₃ℤ)
        [s; inv.(s)]
    end
    half_basis, sizes = Groups.wlmetric_ball(S, radius = half_radius)

    F_sl_3_z = FreeGroup(alphabet(SL₃ℤ))
    e12, e13, e21, e23, e31, e32 = Groups.gens(F_sl_3_z)

    quotient_hom = let source = F_sl_3_z, target = SL₃ℤ
        PropertyT_new.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    relations = [
        e12 * e13 * e12^(-1) * e13^(-1),
        e12 * e32 * e12^(-1) * e32^(-1),
        e13 * e23 * e13^(-1) * e23^(-1),
        e23 * e21 * e23^(-1) * e21^(-1),
        e21 * e31 * e21^(-1) * e31^(-1),
        e31 * e32 * e31^(-1) * e32^(-1),
        e12 * e23 * e12^(-1) * e23^(-1) * e13^(-1),
        e13 * e32 * e13^(-1) * e32^(-1) * e12^(-1),
        e21 * e13 * e21^(-1) * e13^(-1) * e23^(-1),
        e23 * e31 * e23^(-1) * e31^(-1) * e21^(-1),
        e31 * e12 * e31^(-1) * e12^(-1) * e32^(-1),
        e32 * e21 * e32^(-1) * e21^(-1) * e31^(-1),
    ]

    Δ₁x, Iₙ = LowCohomologySOS.spectral_gaps_elements(quotient_hom, relations, half_basis)
    Δ₁x, Iₙ, half_basis
end

Δ₁x_sgap_problem = LowCohomologySOS.sos_problem_matrix(Δ₁x, Iₙ)

warm = nothing

status, warm = PropertyT_new.solve(
    Δ₁x_sgap_problem,
    scs_opt(eps = 1e-9, max_iters = 100, accel=-20),
    warm,
)

@info "Termination status:" status

λ_certified = let (λₐₚ, Qₐₚ) = LowCohomologySOS.get_solution(Δ₁x_sgap_problem)

    @info "Approximated λ: " λₐₚ

    λ_certified = LowCohomologySOS.certify_sos_decomposition(
        Δ₁x,
        Iₙ,
        λₐₚ,
        Qₐₚ,
        half_basis,
        parent(first(Δ₁x))
    )

    @info "Certified λ (interval atithmetic): " λ_certified
    λ_certified
end
