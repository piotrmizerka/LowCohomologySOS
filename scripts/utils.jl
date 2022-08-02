import JuMP
import JuMP.MOI

using Dates
using Serialization
using Logging

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

        if Sys.iswindows()
            log_file = joinpath(logdir, replace("solver_$date.log",":" => "_"))
        else
            log_file = joinpath(logdir, "solver_$date.log")
        end

        status, warm = @time solve(log_file, model, optimizer, warm)

        λ, Q = LowCohomologySOS.get_solution(model)
        solution = Dict(:λ=>λ, :Q=>Q, :warm=>warm)
        
        if Sys.iswindows()
            serialize(joinpath(logdir, replace("solution_$date.sjl", ":" => "_")), solution)
        else
            serialize(joinpath(logdir, "solution.sjl"), solution)
        end

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

##
# Low-level solve

setwarmstart!(model::JuMP.Model, ::Nothing) = model

function setwarmstart!(model::JuMP.Model, warmstart)
    constraint_map = Dict(
        ct => JuMP.all_constraints(model, ct...) for
        ct in JuMP.list_of_constraint_types(model)
    )

    JuMP.set_start_value.(JuMP.all_variables(model), warmstart.primal)

    for (ct, idx) in pairs(constraint_map)
        JuMP.set_start_value.(idx, warmstart.slack[ct])
        JuMP.set_dual_start_value.(idx, warmstart.dual[ct])
    end
    return model
end

function getwarmstart(model::JuMP.Model)
    constraint_map = Dict(
        ct => JuMP.all_constraints(model, ct...) for
        ct in JuMP.list_of_constraint_types(model)
    )

    primal = JuMP.value.(JuMP.all_variables(model))

    slack = Dict(k => JuMP.value.(v) for (k, v) in constraint_map)
    duals = Dict(k => JuMP.dual.(v) for (k, v) in constraint_map)

    return (primal = primal, dual = duals, slack = slack)
end

function solve(m::JuMP.Model, optimizer, warmstart = nothing)

    JuMP.set_optimizer(m, optimizer)
    JuMP.MOIU.attach_optimizer(m)

    m = setwarmstart!(m, warmstart)

    JuMP.optimize!(m)
    Base.Libc.flush_cstdio()

    status = JuMP.termination_status(m)

    return status, getwarmstart(m)
end

function solve(solverlog::String, m::JuMP.Model, optimizer, warmstart = nothing)

    isdir(dirname(solverlog)) || mkpath(dirname(solverlog))

    Base.flush(Base.stdout)
    Base.Libc.flush_cstdio()
    status, warmstart = open(solverlog, "a+") do logfile
        redirect_stdout(logfile) do
            status, warmstart = solve(m, optimizer, warmstart)
            status, warmstart
        end
    end

    return status, warmstart
end
