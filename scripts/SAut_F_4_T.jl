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
const half_radius = 3
sautfN = Groups.SpecialAutomorphismGroup(FreeGroup(N))
RG, S, sizes = @time PropertyT.group_algebra(sautfN, halfradius=half_radius, twisted=true)

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
        StarAlgebras.Basis{UInt32}(@view basis(RG)[1:sizes[half_radius]]),
    )
end

order_unit = Δ = RG(length(S)) - sum(RG(s) for s in S)
eoi = Δ^2 # we don't have to untwist the coeffs since Δ is hermitian

@time sos_problem, P = PropertyT.sos_problem_primal(
    eoi, 
    order_unit,
    wedderburn_decomposition,
    # upper_bound=UPPER_BOUND,
    augmented=true,
    show_progress=true
)

JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-7, max_iters = 100_000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem, P, wedderburn_decomposition)

PropertyT.certify_solution(
    eoi,
    order_unit,
    λ,
    Q,
    halfradius=half_radius,
)
