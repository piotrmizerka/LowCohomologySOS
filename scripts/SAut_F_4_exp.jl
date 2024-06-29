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
relations = LowCohomologySOS.relations(sautfN, F_G, S, true, N, "adj_op")
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
Δ1, I, Adj, Δ1⁻ = LowCohomologySOS.laplacians(sautfN, half_basis_restr, S, sq_adj_op_ = "adj")
Δ1, I, Op, Δ1⁻ = LowCohomologySOS.laplacians(sautfN, half_basis_restr, S, sq_adj_op_ = "op")
RG_star = parent(first(Adj))
Op = LowCohomologySOS.embed.(identity, Op, Ref(RG_star))
I = LowCohomologySOS.embed.(identity, I, Ref(RG_star))
Δ1⁻ = LowCohomologySOS.embed.(identity, Δ1⁻, Ref(RG_star))
sq, adj, op = LowCohomologySOS.sq_adj_op(Δ1⁻, S)
Adj_Op = Adj+adj+100*(Op+op) # KKN had 100 as well :)

# Symmetrize the problem using Wedderburn to accelerate the computations ##############
function wedderburn_data(basis, half_basis, S)
    Σ = PermutationGroups.SymmetricGroup(N)
    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    return constraints_basis, psd_basis, Σ, actions
end
constraints_basis, psd_basis, Σ, action = wedderburn_data(RG_star.basis, half_basis_restr, S)

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Adj_Op, σ, action.alphabet_perm, S) == Adj_Op
    @assert LowCohomologySOS.act_on_matrix(I, σ, action.alphabet_perm, S) == I
end
w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
#################

# Find a numerical spectral gap
sos_problem, P = LowCohomologySOS.sos_problem(Adj_Op, I, w_dec_matrix)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-7, max_iters = 10_000))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
LowCohomologySOS.certify_sos_decomposition(Adj_Op, I, λ, Q, half_basis_restr)
