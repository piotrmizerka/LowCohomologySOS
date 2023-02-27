include(joinpath(@__DIR__, "adj_and_packages.jl"));

constraints_basis, psd_basis, Σ, action = wedderburn_data(basis, half_basis, S);

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Adj, σ, action.alphabet_perm, S) == Adj
    @assert LowCohomologySOS.act_on_matrix(Iₙ, σ, action.alphabet_perm, S) == Iₙ
end

SymbolicWedderburn._int_type(::Type{<:SymbolicWedderburn.InducedActionHomomorphism}) = UInt32

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    @info "Saving Wedderburn decomposition to the file:"
    serialize(joinpath(@__DIR__, "./wedderburn_dec.sjl"), w_dec_matrix)
end
