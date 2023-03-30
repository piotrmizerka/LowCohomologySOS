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
slN = MatrixGroups.SpecialLinearGroup{N}(Int8)
F_sl_N_z = FreeGroup(alphabet(slN))
S = gens(slN)
relations = LowCohomologySOS.relations(slN, F_sl_N_z, S, true, N, "adj")
d₁ = LowCohomologySOS.jacobian_matrix(relations, gens(F_sl_N_z))

# Compute the differentials supported on ball of readius 2 in the standard gen set
S_inv = [S; inv.(S)]
half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2)
RG = LowCohomologySOS.group_ring(slN, half_basis, star_multiplication = false)
quotient_hom = let source = F_sl_N_z, target = slN
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end
d0x = LowCohomologySOS.embed.(Ref(quotient_hom), LowCohomologySOS.d₀(parent(first(d₁)), gens(F_sl_N_z)), Ref(RG))
d1x = LowCohomologySOS.embed.(Ref(quotient_hom), d₁, Ref(RG))

# Restrict attention to the differentials' support
half_basis_restr = [one(slN)]
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
Δ1, I, Δ1⁺, Δ1⁻ = LowCohomologySOS.laplacians(slN, half_basis_restr, S, sq_adj_op_ = "adj")
RG_star = parent(first(Δ1))
I = [i ≠ j ? zero(RG_star) : one(RG_star) for i in 1:length(d0x), j in 1:length(d0x)]
sq, adj, op = LowCohomologySOS.sq_adj_op(Δ1⁻, S)
Adj = Δ1⁺+adj

# Symmetrize the problem using Wedderburn to accelerate the computations ##############
function wedderburn_data(basis, half_basis, S)
    Σ = PermutationGroups.SymmetricGroup(N)
    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    return constraints_basis, psd_basis, Σ, actions
end
constraints_basis, psd_basis, Σ, action = wedderburn_data(parent(first(I)).basis, half_basis_restr, S);
# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    # @assert LowCohomologySOS.act_on_matrix(Δ1, σ, action.alphabet_perm, S) == Δ1
    @assert LowCohomologySOS.act_on_matrix(Adj, σ, action.alphabet_perm, S) == Adj
    @assert LowCohomologySOS.act_on_matrix(I, σ, action.alphabet_perm, S) == I
end
w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
#################

# Find a numerical spectral gap
# sos_problem, P = LowCohomologySOS.sos_problem(Δ1, I, w_dec_matrix)
sos_problem, P = LowCohomologySOS.sos_problem(Adj, I, w_dec_matrix)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-6, max_iters = 30_000))
# JuMP.set_optimizer(sos_problem, cosmo_opt(eps = 1e-7, max_iters = 30_000))
JuMP.optimize!(sos_problem)

# Certify the numerical estimate
λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
# LowCohomologySOS.certify_sos_decomposition(Δ1, I, λ, Q, half_basis_restr)
LowCohomologySOS.certify_sos_decomposition(Adj, I, λ, Q, half_basis_restr)

# Let's try to find somme pattern in the solution
RG = LowCohomologySOS.group_ring(slN,half_basis_restr,star_multiplication = false)
summand_factors = LowCohomologySOS.sos_summand_factors(RG, Q, half_basis_restr, 0.000001, 6);
sos_solution_heuristic = sum(Ref(@interval(1)).*ξ'*ξ for ξ in summand_factors)
Adjx, Ix = LowCohomologySOS.embed.(identity,Adj,Ref(RG)), LowCohomologySOS.embed.(identity,I,Ref(RG))
residual = Adjx-Ref(@interval(λ)).*Ix-sos_solution_heuristic
l1_norm = sum(x -> norm(x, 1), residual)
@interval(λ)-l1_norm
