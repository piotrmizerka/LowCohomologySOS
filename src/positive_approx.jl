function constraints(pm::AbstractMatrix{<:Integer})
    @assert all(i->1≤length(pm), pm)
    cnstrs = [Vector{Int}() for _ in 1:maximum(pm)]
    li = LinearIndices(pm)

    for (idx, k) in pairs(pm)
        push!(cnstrs[k], li[idx])
    end

    Threads.@threads for i in 1:length(cnstrs)
        sort!(cnstrs[i])
    end

    a,b = size(pm)

    return [BinaryMatrix(c, a, b, 1, sorted = true) for c in cnstrs]
end

function sos_problem_matrix(
    M::AbstractMatrix{<:AlgebraElement},
    order_unit::AbstractMatrix{<:AlgebraElement},
    upper_bound::Float64 = Inf,
)
    @assert size(M) == size(order_unit)
    @assert !isempty(M)

    A = parent(first(M))
    @assert all(x -> parent(x) === A, M)
    @assert all(x -> parent(x) === A, order_unit)

    m = size(A.mstructure, 1)
    n = LinearAlgebra.checksquare(M)
    mn = m * n
    result = JuMP.Model()

    JuMP.@variable(result, P[1:mn, 1:mn], Symmetric)
    JuMP.@constraint(result, sdp, P in JuMP.PSDCone())

    JuMP.@variable(result, λ)
    JuMP.@objective(result, Max, λ)

    if upper_bound < Inf
        λ = JuMP.@constraint(result, λ <= upper_bound)
    end

    cnstrs = constraints(A.mstructure)
    @assert length(cnstrs) == length(basis(A))

    for idx in CartesianIndices(M)
        i, j = Tuple(idx)
        Pⁱʲ = @view P[KroneckerDelta{n}(i, j)]

        mij = M[idx]
        uij = order_unit[idx]

        for (A_g, g) in zip(cnstrs, basis(A))
            JuMP.@constraint(result, mij(g) - λ * uij(g) == dot(A_g, Pⁱʲ))
        end
    end
    return result
end

function sos_problem_solution(sos_problem; optimizer)
    JuMP.set_optimizer(sos_problem, optimizer)
    JuMP.optimize!(sos_problem)
    λ, P, termination_status = JuMP.value(sos_problem[:λ]), JuMP.value.(sos_problem[:P]), JuMP.termination_status(sos_problem)

    return λ, P, termination_status
end

# h:Free group --> our group G
function sos_problem_delta_1(
    h,
    relations::AbstractVector{<:FPGroupElement},
    half_basis,
)
    @assert !isempty(relations)
    F = parent(first(relations)) # source of h
    G = parent(h(first(relations))) # target of h

    d₁ = jacobian_matrix(relations)

    RG_ball = group_ring(G, half_basis, star_multiplication = false)

    d₁x = embed.(Ref(h), d₁, Ref(RG_ball))
    d₀x = embed.(Ref(h), d₀(parent(first(d₁)), Groups.gens(F)), Ref(RG_ball))
    Δ₁⁺ = d₁x' * d₁x
    Δ₁⁻ = d₀x * d₀x'
    Δ₁ = Δ₁⁺ + Δ₁⁻

    RG_ball_star = group_ring(G, half_basis, star_multiplication = true)

    Δ₁x = embed.(identity, Δ₁, Ref(RG_ball_star))

    n = length(Groups.gens(F))
    @assert size(Δ₁x, 1) === size(Δ₁x, 2) === n
    Iₙ = [i ≠ j ? zero(RG_ball_star) : one(RG_ball_star) for i in 1:n, j in 1:n]

    Δ₁_sos_problem = sos_problem_matrix(Δ₁x, Iₙ)

    return Δ₁_sos_problem, Δ₁x, Iₙ, half_basis, RG_ball_star
end

# h:Free group --> our group G
function spectral_gaps_approximated(
    h,
    relations::AbstractVector{<:FPGroupElement},
    half_basis;
    optimizer,
)
    Δ₁_sos_problem, Δ₁, Iₙ, half_basisx = sos_problem_delta_1(h, relations, half_basis)
    λ, P, termination_status = sos_problem_solution(Δ₁_sos_problem, optimizer = optimizer)

    return λ, P, termination_status, Δ₁, Iₙ
end
