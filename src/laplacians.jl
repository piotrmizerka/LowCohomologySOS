function laplacians(
    G, # either SL(n,ℤ) or SAut(Fₙ)
    half_basis,
    S; # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    twist_coeffs = true
)
    N = typeof(G) <: MatrixGroups.SpecialLinearGroup ? size(first(gens(G)))[1] : length(G.domain)
    
    @assert (typeof(G) <: MatrixGroups.SpecialLinearGroup && (length(S) == N*(N-1) || length(S) == 2*N*(N-1))) ||
            (typeof(G) <: AutomorphismGroup && (length(S) == 2*N*(N-1) || length(S) == 4*N*(N-1)))

    symmetric_action = (typeof(G) <: MatrixGroups.SpecialLinearGroup && length(S) == N*(N-1)) ||
                       (typeof(G) <: AutomorphismGroup && length(S) == 2*N*(N-1)) ? true : false
    
    if symmetric_action
        F_G = FreeGroup(alphabet(G))
        quotient_hom = let source = F_G, target = G
            Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
        end
    else
        F_G = FreeGroup(length(S))
        quotient_hom = let source = F_G, target = G
            Groups.Homomorphism((i, F, G) -> Groups.word_type(G)(
                [free_group_index(i,F_G, S)]), source, target)
        end
    end

    # check if the quotient homomorphism is defined properly
    @assert length(gens(F_G)) == length(S)
    for i in eachindex(S)
        @assert quotient_hom(gens(F_G,i)) == S[i]
        @assert quotient_hom(gens(F_G,i)^(-1)) == S[i]^(-1)
    end

    relationsx = relations(G, F_G, S, symmetric_action, N)

    return LowCohomologySOS.spectral_gap_elements(quotient_hom, relationsx, half_basis, twist_coeffs = twist_coeffs)
end

# relations for G = SL(n,ℤ), SAut(Fₙ) for symmetric and wreath actions
function relations(
    G,
    F_G::Groups.FreeGroup,
    S, # the generating set for G: either elementary matrices for SL(n,ℤ) or Nielsen transvections for SAut(Fₙ)
    symmetric_action::Bool, # true for the action of Sₙ, false assuming the action of ℤ₂≀Sₙ
    N::Integer
)
    gen_dict = Dict(determine_letter(S[i]) => gens(F_G, i) for i in eachindex(S))

    range_as_list = [i for i in 1:N]
    pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ 1:N
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]

    if typeof(G) <: MatrixGroups.SpecialLinearGroup
        # The presentation taken from the article of Conder et. al.: https://www.jstor.org/stable/2159559#metadata_info_tab_contents 
        e(i,j,ε) = gen_dict[MatrixGroups.ElementaryMatrix{N}(i,j,Int8(ε)*Int8(1))]

        # quadruples_total = [(i,j,k,m) for k ∈ 1:N
        #                         for m ∈ deleteat!(copy(range_as_list), findall(m->m==k,copy(range_as_list)))
        #                         for i ∈ deleteat!(copy(range_as_list), findall(i->i==m,copy(range_as_list))) 
        #                         for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
        # quadruples_wrong_inds = []
        # for ind in eachindex(quadruples_total)
        #     (i,j,k,m) = quadruples_total[ind]
        #     if (i,j) == (k,m)
        #         append!(quadruples_wrong_inds, ind)
        #     end
        # end
        # quadruples_wrong = [quadruples_total[ind] for ind in quadruples_wrong_inds]
        # quadruples = setdiff(quadruples_total, quadruples_wrong)
        
        if symmetric_action
            e(i,j) = e(i,j,1)
            relations = vcat(
                # [e(i,j)*e(k,m)*e(i,j)^(-1)*e(k,m)^(-1) for (i,j,k,m) ∈ quadruples], # for now, let's try to remove these relations
                [e(i,j)*e(i,k)*e(i,j)^(-1)*e(i,k)^(-1) for (i,j,k) ∈ triples],
                [e(i,j)*e(k,j)*e(i,j)^(-1)*e(k,j)^(-1) for (i,j,k) ∈ triples],
                [e(i,j)*e(j,k)*e(i,j)^(-1)*e(j,k)^(-1)*e(i,k)^(-1) for (i,j,k) ∈ triples]
            )
        else # wreath product action
            relations = vcat(
                [e(i,j,ε)*e(i,j,-ε) for (i,j) ∈ pairs for ε ∈ [1,-1]],
                # [e(i,j,ε₁)*e(k,l,ε₂)*e(i,j,-ε₁)*e(k,l,-ε₂) for (i,j,k,l) ∈ quadruples for ε₁ ∈ [1,-1] for ε₂ ∈ [1,-1]],
                [e(i,j,ε₁)*e(i,k,ε₂)*e(i,j,-ε₁)*e(i,k,-ε₂) for (i,j,k) ∈ triples for ε₁ ∈ [1,-1] for ε₂ ∈ [1,-1]],
                [e(i,j,ε₁)*e(k,j,ε₂)*e(i,j,-ε₁)*e(k,j,-ε₂) for (i,j,k) ∈ triples for ε₁ ∈ [1,-1] for ε₂ ∈ [1,-1]],
                [e(i,j,ε₁)*e(j,k,ε₂)*e(i,j,-ε₁)*e(j,k,-ε₂)*e(i,k,-ε₁*ε₂) for (i,j,k) ∈ triples for ε₁ ∈ [1,-1] for ε₂ ∈ [1,-1]]
            )
        end
    else
        # The relations are derived from the Gersten's article, https://www.sciencedirect.com/science/article/pii/0022404984900628,
        # Theorem 2.8. Note that our convention assumes the automorphism composition order reversed with respect to Gersten's.
        # Therefore, the order of letters in the relators had to be reversed as well (earlier, we had to change Gertsen's symbols E_a_b
        # to the ϱ and λ notation used in https://annals.math.princeton.edu/2021/193-2/p03).
        ϱ(i,j, ε) = gen_dict[Groups.Transvection(:ϱ, i, j, ε)]
        λ(i,j, ε) = gen_dict[Groups.Transvection(:λ, i, j, ε)]

        quadruples_1 = [(i,j,k,l) for k ∈ 1:N
                            for l ∈ deleteat!(copy(range_as_list), findall(l->l==k,copy(range_as_list))) 
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i∈[k,l],copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
        quadruples_2 = [(i,j,k,l) for k ∈ 1:N
                            for l ∈ deleteat!(copy(range_as_list), findall(l->l==k,copy(range_as_list)))
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i==l,copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]

        if symmetric_action
            ϱ(i,j) = ϱ(i,j,false)
            λ(i,j) = λ(i,j,false)

            relations = vcat(
                [ϱ(k,l)^(-1)*ϱ(i,j)^(-1)*ϱ(k,l)*ϱ(i,j) for (i,j,k,l) ∈ quadruples_1],
                [λ(k,l)^(-1)*λ(i,j)^(-1)*λ(k,l)*λ(i,j) for (i,j,k,l) ∈ quadruples_1],
                [λ(k,l)^(-1)*ϱ(i,j)^(-1)*λ(k,l)*ϱ(i,j) for (i,j,k,l) ∈ quadruples_2],

                [ϱ(i,k)^(-1)*ϱ(j,k)^(-1)*ϱ(i,j)^(-1)*ϱ(j,k)*ϱ(i,j) for (i,j,k) ∈ triples],
                [ϱ(i,k)*ϱ(j,k)*ϱ(i,j)^(-1)*ϱ(j,k)^(-1)*ϱ(i,j) for (i,j,k) ∈ triples],
                [λ(i,k)^(-1)*λ(j,k)^(-1)*λ(i,j)^(-1)*λ(j,k)*λ(i,j) for (i,j,k) ∈ triples],
                [λ(i,k)*λ(j,k)*λ(i,j)^(-1)*λ(j,k)^(-1)*λ(i,j) for (i,j,k) ∈ triples],
                [ϱ(i,k)*λ(j,k)^(-1)*ϱ(i,j)*λ(j,k)*ϱ(i,j)^(-1) for (i,j,k) ∈ triples],
                [ϱ(i,k)^(-1)*λ(j,k)*ϱ(i,j)*λ(j,k)^(-1)*ϱ(i,j)^(-1) for (i,j,k) ∈ triples],
                [λ(i,k)*ϱ(j,k)^(-1)*λ(i,j)*ϱ(j,k)*λ(i,j)^(-1) for (i,j,k) ∈ triples],
                [λ(i,k)^(-1)*ϱ(j,k)*λ(i,j)*ϱ(j,k)^(-1)*λ(i,j)^(-1) for (i,j,k) ∈ triples]
            )
        else
            relations = vcat(
                [ϱ(i,j,ε)*ϱ(i,j,!ε) for (i,j) ∈ pairs for ε ∈ [true,false]],
                [λ(i,j,ε)*λ(i,j,!ε) for (i,j) ∈ pairs for ε ∈ [true,false]],

                [ϱ(k,l,!ε₂)*ϱ(i,j,!ε₁)*ϱ(k,l,ε₂)*ϱ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_1 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
                [λ(k,l,!ε₂)*λ(i,j,!ε₁)*λ(k,l,ε₂)*λ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_1 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
                [ϱ(k,l,!ε₂)*λ(i,j,!ε₁)*ϱ(k,l,ε₂)*λ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_2 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],
                [λ(k,l,!ε₂)*ϱ(i,j,!ε₁)*λ(k,l,ε₂)*ϱ(i,j,ε₁) for (i,j,k,l) ∈ quadruples_2 for ε₁ ∈ [true,false] for ε₂ ∈ [true,false]],

                [ϱ(i,k)*ϱ(j,k)*ϱ(i,j,true)*ϱ(j,k,true)*ϱ(i,j) for (i,j,k) ∈ triples],
                [ϱ(i,k,true)*ϱ(j,k,true)*ϱ(i,j,true)*ϱ(j,k)*ϱ(i,j) for (i,j,k) ∈ triples],
                [λ(i,k)*λ(j,k)*λ(i,j,true)*λ(j,k,true)*λ(i,j) for (i,j,k) ∈ triples],
                [λ(i,k,true)*λ(j,k,true)*λ(i,j,true)*λ(j,k)*λ(i,j) for (i,j,k) ∈ triples],
                [ϱ(i,k)*λ(j,k,true)*ϱ(i,j)*λ(j,k)*ϱ(i,j,true) for (i,j,k) ∈ triples],
                [ϱ(i,k,true)*λ(j,k)*ϱ(i,j)*λ(j,k,true)*ϱ(i,j,true) for (i,j,k) ∈ triples],
                [λ(i,k)*ϱ(j,k,true)*λ(i,j)*ϱ(j,k)*λ(i,j,true) for (i,j,k) ∈ triples],
                [λ(i,k,true)*ϱ(j,k)*λ(i,j)*ϱ(j,k,true)*λ(i,j,true) for (i,j,k) ∈ triples]
            )
        end
    end

    return relations
end

function determine_letter(g)
    @assert length(word(g)) == 1

    A = alphabet(parent(g))

    return A[first(word(g))]
end

function free_group_index(i::Integer, F_G, S)
    gen_id = (-1, false)
    for j in eachindex(gens(F_G))
        if gens(F_G, j) == F_G([i])
            gen_id = (j, false)
        elseif gens(F_G, j) == inv(F_G([i]))
            gen_id = (j, true)
        end
    end
    
    return gen_id[2] ? word(inv(S[gen_id[1]]))[1] : word(S[gen_id[1]])[1]
end

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