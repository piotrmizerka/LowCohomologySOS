function sln_laplacians(
    slN,
    half_basis,
    S # the generating set for SL(n,ℤ)
)
    @assert typeof(slN) <: MatrixGroups.SpecialLinearGroup

    F_sl_4_z = FreeGroup(alphabet(slN))

    quotient_hom = let source = F_sl_4_z, target = slN
        Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
    end

    elmatrix_gen_dict = Dict(determine_letter(S[i]) => gens(F_sl_4_z, i) for i in eachindex(S))
    N = size(first(gens(slN)))[1]
    e(i,j) = elmatrix_gen_dict[MatrixGroups.ElementaryMatrix{N}(i,j,Int8(1))]

    range_as_list = [i for i in 1:N]
    quadruples_total = [(i,j,k,m) for k ∈ 1:N
                            for m ∈ deleteat!(copy(range_as_list), findall(m->m==k,copy(range_as_list)))
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i==m,copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    quadruples_wrong_inds = []
    for ind in eachindex(quadruples_total)
        (i,j,k,m) = quadruples_total[ind]
        if (i,j) == (k,m)
            append!(quadruples_wrong_inds, ind)
        end
    end
    quadruples_wrong = [quadruples_total[ind] for ind in quadruples_wrong_inds]
    quadruples = setdiff(quadruples_total, quadruples_wrong)
    triples = [(i,j,k) for i ∈ 1:N
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]
    
    # The presentation taken from the article of Conder et. al.: https://www.jstor.org/stable/2159559#metadata_info_tab_contents 
    relations = vcat(
        [e(i,j)*e(k,m)*e(i,j)^(-1)*e(k,m)^(-1) for (i,j,k,m) ∈ quadruples],
        [e(i,j)*e(j,k)*e(i,j)^(-1)*e(j,k)^(-1)*e(i,k)^(-1) for (i,j,k) ∈ triples]
    )

    return LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)
end

function sautfn_laplacians(
    sautfN,
    half_basis,
    S, # the generating set for autfN
    wreath_action # this parameter can be derived from the size of the generating set S (TODO?)
)
    @assert typeof(sautfN) <: AutomorphismGroup

    if !wreath_action
        F_G = FreeGroup(alphabet(sautfN))
        quotient_hom = let source = F_G, target = sautfN
            Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
        end
    else
        F_G = FreeGroup(length(S))
        quotient_hom = let source = F_G, target = sautfN
            Groups.Homomorphism((i, F, G) -> Groups.word_type(G)(
                [free_group_saut_index(i,F_G, S)]), source, target)
        end
    end

    # check if the quotient homomorphism is defined properly
    @assert length(gens(F_G)) == length(S)
    for i in 1:length(S)
        @assert quotient_hom(gens(F_G,i)) == S[i]
    end

    N = length(sautfN.domain)

    transvection_gen_dict = Dict([(determine_letter(S[i]), gens(F_G, i)) for i in 1:length(S)])
    ϱ(i,j, ε) = transvection_gen_dict[Groups.Transvection(:ϱ, i, j, ε)]
    λ(i,j, ε) = transvection_gen_dict[Groups.Transvection(:λ, i, j, ε)]

    range_as_list = [i for i in 1:N]
    quadruples_1 = [(i,j,k,l) for k ∈ 1:N
                            for l ∈ deleteat!(copy(range_as_list), findall(l->l==k,copy(range_as_list))) 
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i∈[k,l],copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    quadruples_2 = [(i,j,k,l) for k ∈ 1:N
                            for l ∈ deleteat!(copy(range_as_list), findall(l->l==k,copy(range_as_list)))
                            for i ∈ deleteat!(copy(range_as_list), findall(i->i==l,copy(range_as_list))) 
                            for j ∈ deleteat!(copy(range_as_list), findall(j->j∈[i,k],copy(range_as_list)))]
    triples = [(i,j,k) for i ∈ 1:N
                    for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list))) 
                    for k ∈ deleteat!(copy(range_as_list), findall(k->k∈[i,j],copy(range_as_list)))]

    # The relations are derived from the Gersten's article, https://www.sciencedirect.com/science/article/pii/0022404984900628,
    # Theorem 2.8. Note that our convention assumes the automorphism composition order reversed with respect to Gersten's.
    # Therefore, the order of letters in the relators had to be reversed as well (earlier, we had to change Gertsen's symbols E_a_b
    # to the ϱ and λ notation used in https://annals.math.princeton.edu/2021/193-2/p03).
    if wreath_action
        pairs = [(i,j) for i ∈ 1:N for j ∈ deleteat!(copy(range_as_list), findall(j->j==i,copy(range_as_list)))]

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
    else
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
    end

    return LowCohomologySOS.spectral_gap_elements(quotient_hom, relations, half_basis)
end

function determine_letter(g)
    @assert length(word(g)) == 1
    
    A = alphabet(parent(g))

    return A[first(word(g))]
end

function free_group_saut_index(i::Integer, F_G, S)
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