using Revise
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "./")))

using LowCohomologySOS
using Groups
using SymbolicWedderburn
using PermutationGroups

include(joinpath(@__DIR__, "./scripts/optimizers.jl"))
include(joinpath(@__DIR__, "./scripts/utils.jl"))

function group_data(half_radius, N, wreath_action)
    SAut_F(n) = Groups.SpecialAutomorphismGroup(FreeGroup(n))
    SAut_F_N = SAut_F(N)

    S_inv = let s = gens(SAut_F_N)
        [s; inv.(s)]
    end
    S = (wreath_action ? S_inv : gens(SAut_F_N))
    basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]
    ℝSAutF_N_star = LowCohomologySOS.group_ring(SAut_F_N, half_basis, star_multiplication = true)

    return SAut_F_N, ℝSAutF_N_star.basis, half_basis, S
end

function wedderburn_data(basis, half_basis, S)
    @time begin
        N = length((parent(first(S))).domain)
        if length(S) == 4*N*(N-1)
            Z_2_wr_S(n) = Groups.Constructions.WreathProduct(PermutationGroups.SymmetricGroup(2), PermutationGroups.SymmetricGroup(n))
            Σ = Z_2_wr_S(N)
        else
            Σ = PermutationGroups.SymmetricGroup(N)
        end
        actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)
        constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    end

    return constraints_basis, psd_basis, Σ, actions
end

const N = 4;
const wreath_action = false;

half_radius = 1;

SAut_F_N, basis, half_basis, S = group_data(half_radius, N, wreath_action)

Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(SAut_F_N, half_basis, S, sq_adj_op_ = "op")
sq, adj, op = LowCohomologySOS.sq_adj_op(Δ₁⁻, S)

Op = Δ₁⁺+op

constraints_basis, psd_basis, Σ, action = wedderburn_data(basis, half_basis, S);

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Op, σ, action.alphabet_perm, S) == Op
    @assert LowCohomologySOS.act_on_matrix(Iₙ, σ, action.alphabet_perm, S) == Iₙ
end

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

sos_problem, P = LowCohomologySOS.sos_problem(Op, Iₙ, w_dec_matrix, 0.05)
JuMP.set_optimizer(sos_problem, scs_opt(eps = 1e-9, max_iters = 5000))
JuMP.optimize!(sos_problem)

λ, Q = LowCohomologySOS.get_solution(sos_problem, P, w_dec_matrix)
LowCohomologySOS.certify_sos_decomposition(Op, Iₙ, λ, Q, half_basis)