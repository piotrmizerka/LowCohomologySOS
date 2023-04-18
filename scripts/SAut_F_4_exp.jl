using Revise
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)
using LowCohomologySOS
using Groups
using JuMP
using SparseArrays
using PermutationGroups
using SymbolicWedderburn
using IntervalArithmetic
include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

# Compute Jacobian in the free group ring
const N = 4
sautfN = Groups.SpecialAutomorphismGroup(FreeGroup(N))
F_G = FreeGroup(alphabet(sautfN))
S = gens(sautfN)
relations = LowCohomologySOS.relations(sautfN, F_G, S, true, N, "adj")
d₁ = LowCohomologySOS.jacobian_matrix(relations, gens(F_G))

# Compute the differentials supported on ball of readius 2 in the standard gen set
S_inv = [S; inv.(S)]
half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2)
RG = LowCohomologySOS.group_ring(sautfN, half_basis, star_multiplication = false)
quotient_hom = let source = F_G, target = sautfN
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end
d0x = LowCohomologySOS.embed.(Ref(quotient_hom), LowCohomologySOS.d₀(parent(first(d₁)), gens(F_G)), Ref(RG))
d1x = LowCohomologySOS.embed.(Ref(quotient_hom), d₁, Ref(RG))

# Restrict attention to the differentials' support
half_basis_restr = [one(sautfN)]
for j in eachindex(d0x)
    for i in SparseArrays.nonzeroinds(d0x[j].coeffs)
        push!(half_basis_restr, RG.basis[i])
    end
end
for j in eachindex(d1x)
    for i in SparseArrays.nonzeroinds(d1x[j].coeffs)
        push!(half_basis_restr, RG.basis[i])
    end
end
half_basis_restr = unique!([half_basis_restr; inv.(half_basis_restr)])

# Compute Laplacians over RG over the differentials' support
# Δ1, I, Δ1⁺, Δ1⁻ = LowCohomologySOS.laplacians(sautfN, half_basis_restr, S, sq_adj_op_ = "adj")
Δ1, I, Δ1⁺, Δ1⁻ = LowCohomologySOS.laplacians(sautfN, half_basis_restr, S, sq_adj_op_ = "all")
sq, adj, op = LowCohomologySOS.sq_adj_op(Δ1⁻, S)
Adj = Δ1⁺+adj

# Try further diminishing of the support and adding ball3 elements given by M. Nitsche,
# see https://arxiv.org/abs/2009.05134, https://zenodo.org/record/7065232#.ZD5dm-zMK3I, and https://github.com/MartinNitsche/AutF4-Property-T
function small_support(d1, d2, Δ; cutoff = 0.08)
    @assert parent(first(d1)) == parent(first(d2)) == parent(first(Δ))
    RG = parent(first(Δ))
    G = parent(first(RG.basis))
    half_basis_restr = [one(G)]
    for j in eachindex(d1)
        for i in SparseArrays.nonzeroinds(d1[j].coeffs)
            push!(half_basis_restr, RG.basis[i])
        end
    end
    for j in eachindex(d2)
        for i in SparseArrays.nonzeroinds(d2[j].coeffs)
            push!(half_basis_restr, RG.basis[i])
        end
    end
    half_basis_restr = unique!([half_basis_restr; inv.(half_basis_restr)])
    delta_support = Dict([h^(-1)*g => false for g in half_basis_restr for h in half_basis_restr])
    for entry in Δ
        for i in SparseArrays.nonzeroinds(entry.coeffs)
            delta_support[RG.basis[i]] = true
        end
    end
    temp_result = [one(G)]
    for g in half_basis_restr
        counter = count([delta_support[h^(-1)*g] for h in half_basis_restr])
        if counter/length(half_basis_restr) >= cutoff
            push!(temp_result,g)
        end
    end
    Σ = PermutationGroups.SymmetricGroup(length(G.domain))
    A = alphabet(G)
    alphabet_perm = LowCohomologySOS.AlphabetPermutation(
        Dict(
            σ => PermutationGroups.Perm([A[LowCohomologySOS._conj(l, σ)] for l in A.letters]) for
            σ in Σ
        ),
    )
    result = [one(G)]
    for g in temp_result
        for σ in Σ
            push!(result,SymbolicWedderburn.action(alphabet_perm,σ,g))
        end
    end
    return unique!([result; inv.(result)])
end
function transvection_from_string(gen_str)
    splitted = split(gen_str,"/")
    i, j = parse(Int,splitted[1]), parse(Int,splitted[2])
    if i>0 && j>0
        return (Groups.Transvection(:ϱ, i, j, false), false)
    elseif i<0 && j<0
        return (Groups.Transvection(:λ, -i, -j, false), false)
    elseif i<0 && j>0
        return (Groups.Transvection(:λ, -i, j, false), true)
    else
        return (Groups.Transvection(:ϱ, i, -j, false), true)
    end
end
function ball_3_elts(G, path_3_ball)
    read = open(path_3_ball)
    lines = readlines(read)
    trans_dict = Dict(LowCohomologySOS.determine_letter(s) => s for s in S)
    temp_result = [one(G)]
    for line in lines
        transvections = split(line,", ")
        pop!(transvections)
        if length(transvections) == 3 && count("-",transvections[3])%2==0
            gen1_str, gen2_str = transvections[1], transvections[2]
            trans1, inv1 = transvection_from_string(gen1_str)
            trans2, inv2 = transvection_from_string(gen2_str)
            gen1 = (inv1 ? trans_dict[trans1]^(-1) : trans_dict[trans1])
            gen2 = (inv2 ? trans_dict[trans2]^(-1) : trans_dict[trans2])
            push!(temp_result,gen1*gen2)
        end
    end
    close(read)
    Σ = PermutationGroups.SymmetricGroup(length(G.domain))
    A = alphabet(G)
    alphabet_perm = LowCohomologySOS.AlphabetPermutation(
        Dict(
            σ => PermutationGroups.Perm([A[LowCohomologySOS._conj(l, σ)] for l in A.letters]) for
            σ in Σ
        ),
    )
    result = [one(G)]
    for g in temp_result
        for σ in Σ
            push!(result,SymbolicWedderburn.action(alphabet_perm,σ,g))
        end
    end
    return unique!([result; inv.(result)])
end
additional_support = ball_3_elts(sautfN,"/Users/piotrmizerka/Desktop/postdoc_warsaw/support.dat") # the path to the downloaded (from Zenodo) 3-ball support of M. Nitsche, has to be changed
RG_restr = LowCohomologySOS.group_ring(sautfN, half_basis_restr, star_multiplication = false)
d0 = LowCohomologySOS.embed.(identity,d0x,Ref(RG_restr));
d1 = LowCohomologySOS.embed.(identity,d1x,Ref(RG_restr));
Δ1 = LowCohomologySOS.embed.(identity,Δ1,Ref(RG_restr));
support_restr = small_support(d0,d1,Δ1,cutoff = 0.05)
support = unique!([support_restr;additional_support])
RG_star = LowCohomologySOS.group_ring(sautfN,support,star_multiplication = true)
RG_star.basis
Δ1 = LowCohomologySOS.embed.(identity,Δ1,Ref(RG_star));
I = LowCohomologySOS.embed.(identity,I,Ref(RG_star));

# Symmetrize the problem using Wedderburn to accelerate the computations ##############
function wedderburn_data(basis, half_basis, S)
    Σ = PermutationGroups.SymmetricGroup(N)
    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    return constraints_basis, psd_basis, Σ, actions
end
constraints_basis, psd_basis, Σ, action = wedderburn_data(RG_star.basis, support, S)
# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    # @assert LowCohomologySOS.act_on_matrix(Adj, σ, action.alphabet_perm, S) == Adj
    @assert LowCohomologySOS.act_on_matrix(Δ1, σ, action.alphabet_perm, S) == Δ1
    @assert LowCohomologySOS.act_on_matrix(I, σ, action.alphabet_perm, S) == I
end
w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
#################

# Find a numerical spectral gap
# sos_problem, P = LowCohomologySOS.sos_problem(Adj, I, w_dec_matrix)
sos_problem, P = LowCohomologySOS.sos_problem(Δ1, I, w_dec_matrix)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-5, max_iters = 100_000))
# JuMP.set_optimizer(sos_problem, cosmo_opt(eps = 1e-7, max_iters = 30_000))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
# LowCohomologySOS.certify_sos_decomposition(Adj, I, λ, Q, support_restr)
LowCohomologySOS.certify_sos_decomposition(Δ1, I, λ, Q, support_restr)
