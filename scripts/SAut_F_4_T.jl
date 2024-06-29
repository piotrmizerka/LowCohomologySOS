# The code inspired with the code from here: https://github.com/kalmarek/PropertyT.jl/blob/0127d05594b83dd03b2908ade48368cf626e8620/scripts/SpN_Adj.jl#L67

using Revise
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using Groups
using JuMP
using LowCohomologySOS
using PermutationGroups
using PropertyT
using StarAlgebras
using SymbolicWedderburn
include(joinpath(@__DIR__, "optimizers.jl"))

const N = 4
const half_radius = 2
sautfN = Groups.SpecialAutomorphismGroup(FreeGroup(N))
S = [gens(sautfN);inv.(gens(sautfN))]
ball2, sizes = Groups.wlmetric_ball(S,radius=half_radius) # actually, no need to compute ball2 here,
# since it turns out to be contained in the additional_support from M. Nitsche
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
        elt = one(G)
        for transx in transvections
            trans, inv = transvection_from_string(transx)
            gen = (inv ? trans_dict[trans]^(-1) : trans_dict[trans])
            elt *= gen
        end
        push!(temp_result,elt)
    end
    close(read)
    P = PermGroup(perm"(1,2)", Perm(circshift(1:N, -1)))
    Σ = Groups.Constructions.WreathProduct(PermGroup(perm"(1,2)"), P)
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
additional_support = ball_3_elts(sautfN,joinpath(@__DIR__, "./M_Nitsche_support.dat")) # the path to the downloaded (from Zenodo) 3-ball support of M. Nitsche, has to be changed
support = unique!([ball2;additional_support])
RG = LowCohomologySOS.group_ring(sautfN,support,star_multiplication=true)

wedderburn_decomposition = let RG = RG, N = N
    G = StarAlgebras.object(RG)
    P = PermGroup(perm"(1,2)", Perm(circshift(1:N, -1)))
    Σ = Groups.Constructions.WreathProduct(PermGroup(perm"(1,2)"), P)
    act = PropertyT.action_by_conjugation(G, Σ)
    @info "Computing WedderburnDecomposition"
    wdfl = @time SymbolicWedderburn.WedderburnDecomposition(
        Float64,
        Σ,
        act,
        basis(RG),
        StarAlgebras.Basis{UInt32}(@view basis(RG)[1:length(support)]),
    )
end

order_unit = Δ = RG(length(S)) - sum(RG(s) for s in S)
eoi = Δ^2 # we don't have to untwist the coeffs since Δ is hermitian

@time sos_problem, P = PropertyT.sos_problem_primal(
    eoi, 
    order_unit,
    wedderburn_decomposition,
    upper_bound=0.1,
    augmented=true,
    show_progress=true
)

JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-7, max_iters = 1_000_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem, P, wedderburn_decomposition)

PropertyT.certify_solution(
    eoi,
    order_unit,
    λ,
    Q,
    halfradius=half_radius,
)
