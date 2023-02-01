using Revise
using LinearAlgebra
ENV["JULIA_NUM_THREADS"] = Sys.CPU_THREADS÷2
LinearAlgebra.BLAS.set_num_threads(Sys.CPU_THREADS÷2)

using LowCohomologySOS
using Groups
using SymbolicWedderburn
using PermutationGroups

include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

function group_data(half_radius, N, wreath_action)
    SL(n, R) = MatrixGroups.SpecialLinearGroup{n}(R)
    slN = SL(N, Int8)

    S_inv = let s = gens(slN)
        [s; inv.(s)]
    end
    S = (wreath_action ? S_inv : gens(slN))
    basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2*half_radius)
    half_basis = basis[1:sizes[half_radius]]
    RslN_star = LowCohomologySOS.group_ring(slN, half_basis, star_multiplication = true)

    return slN, RslN_star.basis, half_basis, S
end

function wedderburn_data(basis, half_basis, S)
    @time begin
        N = size(first(S))[1]
        if length(S) == 2*N*(N-1)
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
const N = 3;
const wreath_action = false;

slN, basis, half_basis, S = group_data(half_radius, N, wreath_action)

Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.laplacians(slN, half_basis, S)

constraints_basis, psd_basis, Σ, action = wedderburn_data(basis, half_basis, S);

# there is no point of finding a solution if we don't provide invariant matrix
for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Δ₁, σ, action.alphabet_perm, S) == Δ₁
    @assert LowCohomologySOS.act_on_matrix(Iₙ, σ, action.alphabet_perm, S) == Iₙ
end

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    Δ₁_sos_problem = LowCohomologySOS.sos_problem(
        Δ₁, 
        Iₙ,
        w_dec_matrix,
        0.7
    )
end

slN_data = (
    M = Δ₁,
    order_unit = Iₙ,
    half_basis = half_basis,
    RG = parent(first(Δ₁)),
)

solve_in_loop(
    Δ₁_sos_problem,
    w_dec_matrix,
    logdir = "./logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 20_000),
    data = slN_data
)

