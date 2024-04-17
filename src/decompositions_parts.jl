function sq_adj_op(
    Δ₁⁻,
    S # generating set indexing Δ₁⁻
)
    RG = parent(first(Δ₁⁻))
    sln = parent(first(RG.basis))
    sq_pairs = []
    adj_pairs = []
    op_pairs = []
    A = alphabet(sln)
    for s in eachindex(S)
        for t in eachindex(S)
            s_i, s_j = A[word(S[s])[1]].i, A[word(S[s])[1]].j
            t_i, t_j = A[word(S[t])[1]].i, A[word(S[t])[1]].j
            if length(intersect!([s_i,s_j],[t_i,t_j])) == 2
                push!(sq_pairs,(s,t))
            elseif length(intersect!([s_i,s_j],[t_i,t_j])) == 1
                push!(adj_pairs,(s,t))
            else
                push!(op_pairs,(s,t))
            end
        end
    end
    sq = [(i,j) in sq_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    adj = [(i,j) in adj_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]
    op = [(i,j) in op_pairs ? Δ₁⁻[i,j] : zero(RG) for i in eachindex(S), j in eachindex(S)]

    @assert sq+adj+op == Δ₁⁻

    return sq, adj, op
end