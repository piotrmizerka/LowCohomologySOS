include(joinpath(@__DIR__, "adj_and_packages.jl"));

@time begin
    @info "Reading Wedderburn decomposition from file:"
    w_dec_matrix = deserialize(joinpath(@__DIR__, "./wedderburn_dec.sjl"))
end

@time begin
    @info "Defining the SOS problem:"
    sos_problem = LowCohomologySOS.sos_problem(
        Adj, 
        Iₙ,
        w_dec_matrix,
        0.7
    )
end

@time begin
    @info "Saving SOS problem to the file:"
    serialize(joinpath(@__DIR__, "./sos_problem.sjl"), sos_problem)
end
