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

function wedderburn_data(basis, half_basis, S, N, wreath_action)
    @time begin
        if wreath_action
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

function determine_transvection(g)
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

const half_radius = 2;
const N = 4;
const wreath_action = false;

SAut_F_N, basis, half_basis, S = group_data(half_radius, N, wreath_action)

Δ₁, Iₙ, Δ₁⁺, Δ₁⁻ = LowCohomologySOS.sautfn_laplacians(SAut_F_N, half_basis, S, wreath_action)

constraints_basis, psd_basis, Σ, action = wedderburn_data(basis, half_basis, S, N, wreath_action);

@time "\tWedderburn total" begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    Δ₁_sos_problem = LowCohomologySOS.sos_problem(
        Δ₁, 
        Iₙ,
        w_dec_matrix,
        length(collect(Σ)),
        1.0
    )
end

SAut_F_N_data = (
    M = Δ₁,
    order_unit = Iₙ,
    half_basis = half_basis,
    RG = parent(first(Δ₁)),
)

solve_in_loop(
    Δ₁_sos_problem,
    w_dec_matrix,
    logdir = "./LowCohomologySOS/logs",
    optimizer = scs_opt(eps = 1e-9, max_iters = 10_000),
    data = SAut_F_N_data
)
