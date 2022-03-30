using LowCohomologySOS
using Groups
using PropertyT_new
using JuMP
using SCS

function scs_opt(;
    accel = 10,
    alpha = 1.5,
    eps = 1e-9,
    max_iters = 10_000,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        SCS.Optimizer,
        "acceleration_lookback" => accel,
        "acceleration_interval" => max(abs(accel), 1),
        "alpha" => alpha,
        "eps_abs" => eps,
        "eps_rel" => eps,
        "linear_solver" => SCS.DirectSolver,
        "max_iters" => max_iters,
        "warm_start" => true,
        "verbose" => verbose,
    )
end

SL₃ℤ_spectral_gaps = let half_radius = 1
    SL(n, R) = PropertyT_new.SpecialLinearGroup{n}(R)
    SL₃ℤ = SL(3, Int8)

    λ, P, termination_status = LowCohomologySOS.property_t_conjugated_approx(SL₃ℤ, half_radius, optimizer = scs_opt(eps = 1e-5, max_iters = 100_000))
end
