# We want to find alpha and beta such that
# αd₁^*d₁+βd₀d₀^* = symmetrized induced Laplacians from SL(3,Z) to SL(4,Z)

using StarAlgebras
using JuMP

function alpha_beta_problem(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    M::AbstractMatrix{<:AlgebraElement},
    alpha_ge_beta = true
)
    @assert size(A) == size(B) == size(M)
    @assert !isempty(A) && !isempty(B) && !isempty(M)

    RG = parent(first(M))

    @assert all(x -> parent(x) === RG, A)
    @assert all(x -> parent(x) === RG, B)
    @assert all(x -> parent(x) === RG, M)

    result = JuMP.Model()

    JuMP.@variable(result, α)
    JuMP.@variable(result, β)

    # The perfect situation would be in the case α = β
    if alpha_ge_beta 
        JuMP.@constraint(result, α >= β)
        JuMP.@objective(result, Max, β-α)
    else
        JuMP.@constraint(result, β >= α)
        JuMP.@objective(result, Max, α-β)
    end

    for idx in CartesianIndices(M)
        aij = A[idx]
        bij = B[idx]
        mij = M[idx]

        for g in basis(RG)
            JuMP.@constraint(result, α *aij(g)+β*bij(g) == mij(g))
        end
    end

    return result
end

include(joinpath(@__DIR__, "optimizers.jl"))

function alpha_beta_adjust(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    M::AbstractMatrix{<:AlgebraElement},
    alpha_ge_beta = true,
    optimizer = scs_opt(eps = 1e-9, max_iters = 10_000)
)
    problem = alpha_beta_problem(A, B, M, alpha_ge_beta)

    JuMP.set_optimizer(problem, optimizer)
    JuMP.optimize!(problem)

    return JuMP.value.(problem[:α]), JuMP.value.(problem[:β])
end

using Groups

function cyclic_group(n::Integer)
    A = Alphabet([:a, :A], [2, 1])
    F = FreeGroup(A)
    a, = Groups.gens(F)
    e = one(F)
    Cₙ = FPGroup(F, [a^n => e])

    return Cₙ
end

G = cyclic_group(3)

using LowCohomologySOS

RG =  LowCohomologySOS.group_ring(G, 1)

a, = Groups.gens(G)
A = [one(RG) RG(a); one(RG) one(RG)]
B = [one(RG) zero(RG); zero(RG) one(RG)]
M = [2*one(RG) RG(a); one(RG) 2*one(RG)]

alpha_beta_adjustx = alpha_beta_adjust(A, B, M)

alpha_beta_adjustx = alpha_beta_adjust(A, B, M, false)