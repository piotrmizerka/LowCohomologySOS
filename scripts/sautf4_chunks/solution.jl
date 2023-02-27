include(joinpath(@__DIR__, "adj_and_packages.jl"));

@time begin
    @info "Reading SOS problem from the file:"
    sos_problem = deserialize(joinpath(@__DIR__, "./sos_problem.sjl"))
end

@time begin
    @info "Reading Wedderburn decomposition from file:"
    w_dec_matrix = deserialize(joinpath(@__DIR__, "./wedderburn_dec.sjl"))
end

SAut_F_N_data = (
    M = Adj,
    order_unit = Iₙ,
    half_basis = half_basis
)

solve_in_loop(
    sos_problem,
    w_dec_matrix,
    logdir = "./logs_sautf4_adj",
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = SAut_F_N_data
)
