using SCS
using COSMO

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

function cosmo_opt(;
    accel = 10,
    eps = 1e-9,
    max_iters = 10_000,
    verbose = true,
)
    return JuMP.optimizer_with_attributes(
        COSMO.Optimizer,
        "accelerator" => with_options(AndersonAccelerator, mem = accel),
        "eps_abs" => eps,
        "eps_rel" => eps^2, 
        "max_iter" => max_iters, 
        "verbose" => verbose,
    )
end
