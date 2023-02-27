include(joinpath(@__DIR__, "adj_and_packages.jl"));

solution = deserialize(joinpath(@__DIR__, "solution.sjl"))

@time begin
    @info "Certifying the solution loaded from the file:"
    result = LowCohomologySOS.certify_sos_decomposition(
        Adj,
        Iₙ,
        solution[:λ],
        solution[:Q],
        half_basis
    )
end

if result[1]
    @info "Adj part certified with λ = "*string(result[2].lo)
end