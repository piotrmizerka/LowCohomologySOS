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
include(joinpath(@__DIR__, "optimizers.jl"))
include(joinpath(@__DIR__, "utils.jl"))

const N = 4

sautfN = Groups.SpecialAutomorphismGroup(FreeGroup(N))
F_G = FreeGroup(alphabet(sautfN))
S = gens(sautfN)

# relations_op = LowCohomologySOS.relations(sautfN, F_G, S, true, N, "op")
# relations_sq = LowCohomologySOS.relations(sautfN, F_G, S, true, N, "sq")
# relations = union(relations_op,relations_sq)
relations = LowCohomologySOS.relations(sautfN, F_G, S, true, N, "all")

d₁ = LowCohomologySOS.jacobian_matrix(relations, gens(F_G))

S_inv = let s = gens(sautfN)
    [s; inv.(s)]
end
half_basis, sizes = Groups.wlmetric_ball(S_inv, radius = 2)
RG = LowCohomologySOS.group_ring(sautfN, half_basis, star_multiplication = false)
quotient_hom = let source = F_G, target = sautfN
    Groups.Homomorphism((i, F, G) -> Groups.word_type(G)([i]), source, target)
end
d0x = LowCohomologySOS.embed.(Ref(quotient_hom), LowCohomologySOS.d₀(parent(first(d₁)), gens(F_G)), Ref(RG))
d1x = LowCohomologySOS.embed.(Ref(quotient_hom), d₁, Ref(RG))
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
RGx = LowCohomologySOS.group_ring(sautfN, half_basis_restr, star_multiplication = false)
d0xx = LowCohomologySOS.embed.(identity, d0x, Ref(RGx))
d1xx = LowCohomologySOS.embed.(identity, d1x, Ref(RGx))
Δ1⁺x = d1xx' * d1xx
Δ1⁻x = d0xx * d0xx'
RG_star = LowCohomologySOS.group_ring(sautfN, half_basis_restr, star_multiplication = true)

Δ1⁺ = LowCohomologySOS.embed.(identity, Δ1⁺x, Ref(RG_star))
Δ1⁻ = LowCohomologySOS.embed.(identity, Δ1⁻x, Ref(RG_star))
Δ1 = Δ1⁺+Δ1⁻
I = [i ≠ j ? zero(parent(first(Δ1⁺))) : one(parent(first(Δ1⁺))) for i in 1:length(d0x), j in 1:length(d0x)]

function wedderburn_data(basis, half_basis, S)
    Σ = PermutationGroups.SymmetricGroup(N)
    actions = LowCohomologySOS.WedderburnActions(alphabet(parent(first(S))), Σ, LowCohomologySOS._conj, S, basis)
    constraints_basis, psd_basis = LowCohomologySOS.matrix_bases(basis, half_basis, S)
    return constraints_basis, psd_basis, Σ, actions
end

using PermutationGroups

constraints_basis, psd_basis, Σ, action = wedderburn_data(RG_star.basis, half_basis_restr, S)

for σ in Σ
    @assert LowCohomologySOS.act_on_matrix(Δ1, σ, action.alphabet_perm, S) == Δ1
    @assert LowCohomologySOS.act_on_matrix(I, σ, action.alphabet_perm, S) == I
end

using SymbolicWedderburn

@time begin
    @info "Wedderburn:"
    w_dec_matrix = SymbolicWedderburn.WedderburnDecomposition(Float64, Σ, action, constraints_basis, psd_basis)
end

@time begin
    sos_problem = LowCohomologySOS.sos_problem(
        Δ1, 
        I,
        w_dec_matrix,
        0.05
    )
end

JuMP.set_optimizer(sos_problem[1], scs_opt(eps = 1e-9, max_iters = 20_000))
JuMP.optimize!(sos_problem[1])
λ, Q = LowCohomologySOS.get_solution(sos_problem[1],sos_problem[2],w_dec_matrix)
LowCohomologySOS.certify_sos_decomposition(Δ1, I, λ, Q, half_basis_restr)