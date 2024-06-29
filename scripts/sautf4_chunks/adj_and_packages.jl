using Revise
using Pkg
Pkg.activate(normpath(joinpath(@__DIR__, "../../")))
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
using SymbolicWedderburn
using PermutationGroups

include(joinpath(@__DIR__, "../optimizers.jl"))
include(joinpath(@__DIR__, "../utils.jl"))

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

const half_radius = 2;
const N = 4;
const wreath_action = true;

SAut_F_N, basis, half_basis, S = group_data(half_radius, N, wreath_action)

Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(SAut_F_N, half_basis, S, sq_adj_op_ = "adj")
sq, adj, op = LowCohomologySOS.sq_adj_op(Δ₁⁻, S)

Adj = Δ₁⁺+adj