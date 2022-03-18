using LowCohomologySOS
using Groups
using PropertyT_new
using JuMP
using SCS
using Revise

function scs_opt(;
    eps = 1e-5,
    acceleration_lookback = 0,
    max_iters = 20_000,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "eps" => eps,
        "acceleration_lookback" => acceleration_lookback,
        "max_iters" => max_iters,
        "verbose" => verbose,
    )
end

SL₃ℤ_spectral_gaps = let half_radius = 1
    SL(n, R) = PropertyT_new.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    λ, P, termination_status = LowCohomologySOS.property_t_conjugated_approx(SL₃ℤ, half_radius, optimizer = scs_opt(eps = 1e-5, max_iters = 100_000))
end
