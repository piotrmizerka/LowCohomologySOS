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

using Dates
using Serialization
using Logging

# This is the analogue of "solve_in_loop" function from PropertyT_new package of M. Kaluba.
# This version is customized for the matrix case.
function solve_in_loop(model; logdir, optimizer, data)
    @info "logging to $logdir"
    status = MOI.UNKNOWN_RESULT_STATUS
    warm = try
        solution = deserialize(joinpath(logdir, "solution.sjl"))
        warm = solution[:warm]
        @info "trying to warm-start model..."
        warm
    catch
        nothing
    end
    old_lambda = 0.0
    while status != MOI.OPTIMAL
        date = now()
        log_file = joinpath(logdir, "solver_$date.log")
        status, warm = @time PropertyT_new.solve(log_file, model, optimizer, warm)

        λ = JuMP.value(model[:λ])
        Q = real.(sqrt(JuMP.value.(model[:P])))
        solution = Dict(:λ=>λ, :Q=>Q, :warm=>warm)
        serialize(joinpath(logdir, "solution_$date.sjl"), solution)
        serialize(joinpath(logdir, "solution.sjl"), solution)

        flag, certified_λ = open(log_file, append=true) do io
            with_logger(SimpleLogger(io)) do
                LowCohomologySOS.certify_sos_decomposition(
                    data.M, 
                    data.order_unit, 
                    λ, 
                    Q, 
                    data.half_basis, 
                )
            end
        end

        if flag == true && certified_λ.lo ≥ 0
            @info "Certification done with λ = $certified_λ"
            return certified_λ.lo
        end

        rel_change = abs(certified_λ.lo - old_lambda)/(abs(certified_λ.lo) + abs(old_lambda))
        @info "λ = $λ" certified_λ.lo, rel_change

        old_lambda = certified_λ.lo

        if rel_change < 1e-9
            @info "No progress detected, breaking"
            break
        end
    end

    return status == MOI.OPTIMAL ? old_lambda : NaN
end