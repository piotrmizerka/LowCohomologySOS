# The function below concerns definition of constraints defined for the semi-positive definite matrix P to be computed.
# For each entry value, it creates a vector of linear indices of the matrix pm on which this value occurs.
# The function returns the vector of vectors as above indexed by pm matrix's values.
function constraints(pm::AbstractMatrix{<:Integer})
    cnstrs = [Vector{Int}() for _ in 1:maximum(pm)]
    li = LinearIndices(CartesianIndices(size(pm)))

    for i in eachindex(pm)
        push!(cnstrs[pm[i]], li[i])
    end

    return cnstrs
end

# As the constraints function, this function also concerns defining the constraints arising from matrix P.
# More precisely, it defines constraints arising from the (row_id,column_id)-entry of the matrix equation for our problem.
# We can apply the constraints function written by M. Kaluba to define the constraints arising from each entry.
# Order of linear indices for matrices which has to be applied: column order (from top to bottom and from Seattle to Miami).
function entry_constraint(
    cnstrs,
    row_id,
    column_id,
    constraind_id,
    half_radius,
    generators_number,
)
    B =
        (column_id - 1) * half_radius^2 * generators_number +
        (row_id - 1) * half_radius
    result = copy(cnstrs[constraind_id])
    for l in 1:length(cnstrs[constraind_id])
        summand =
            (cnstrs[constraind_id][l] % half_radius != 0) ?
            cnstrs[constraind_id][l] % half_radius : half_radius
        factor = cnstrs[constraind_id][l] - summand
        result[l] = B + factor * generators_number + summand
    end

    return result
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
        mij = M[idx]
        u = StarAlgebras.coeffs(order_unit[idx])
        JuMP.@constraint(
            result,
            [k = 1:length(cnstrs)],
            mij[k] - λ * u[k] ==
            sum(P[p] for p in entry_constraint(cnstrs, Tuple(idx)..., k, m, n))
        )
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
function spectral_gaps_approximated(
    h,
    relations::AbstractVector{<:FPGroupElement},
    half_basis;
    optimizer,
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
    λ, P, termination_status = sos_problem_solution(Δ₁_sos_problem, optimizer = optimizer)

    return λ, P, termination_status, RG_ball_star, Δ₁x, Iₙ
end
