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
# Order of linear indices for matrices which has to be applied: column snake (from Seattle to Miami).
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
    M::AbstractMatrix,
    order_unit,
    upper_bound::Float64 = Inf,
)
    @assert !isempty(M)
    underlying_group_ring = parent(first(M))
    m = size(underlying_group_ring.mstructure, 1)
    n = size(M, 1)
    mn = m * n
    result = JuMP.Model()

    JuMP.@variable(result, P[1:mn, 1:mn], Symmetric)
    JuMP.@constraint(result, sdp, P in JuMP.PSDCone())

    if upper_bound < Inf
        λ = JuMP.@variable(result, λ <= upper_bound)
    else
        λ = JuMP.@variable(result, λ)
    end

    cnstrs = constraints(underlying_group_ring.mstructure)

    for i in 1:n
        for j in 1:n
            mij = StarAlgebras.coeffs(M[i, j])
            u = StarAlgebras.coeffs(order_unit[i, j])
            @assert length(cnstrs) == length(mij) == length(u)
            JuMP.@constraint(
                result,
                [k = 1:length(cnstrs)],
                mij[k] - λ * u[k] ==
                sum(P[entry_constraint(cnstrs, i, j, k, m, n)])
            )
        end
    end

    JuMP.@objective(result, Max, λ)

    return result
end

function sos_problem_solution(sos_problem; optimizer)
    JuMP.set_optimizer(sos_problem, optimizer)
    JuMP.optimize!(sos_problem)
    λ, P = JuMP.value(sos_problem[:λ]), JuMP.value.(sos_problem[:P])

    return λ, P
end

# h:Free group --> our group G
function spectral_gaps_approximated(
    h::Function,
    relations,
    half_basis;
    optimizer,
)
    @assert !isempty(relations)
    F = parent(first(relations))
    G = parent(h(first(relations)))

    d₁ = jacobian_matrix(relations)
    d₀x = d₀(parent(rand(d₁)), Groups.gens(F))
    Δ₁⁺ = d₁' * d₁
    Δ₁⁻ = d₀x * d₀x'
    Δ₁ = Δ₁⁺ + Δ₁⁻

    RG_ball_star = group_ring(G, half_basis, true)

    Δ₁x = embed.(h, Δ₁, Ref(RG_ball_star))

    n = length(Groups.gens(F))
    @assert size(Δ₁x, 1) === size(Δ₁x, 2) === n
    Iₙ = reshape([zero(RG_ball_star) for i in 1:(n*n)], n, n)
    for i in 1:n
        Iₙ[i, i] = one(RG_ball_star)
    end

    Δ₁_sos_problem = sos_problem_matrix(Δ₁x, Iₙ)
    λ, P = sos_problem_solution(Δ₁_sos_problem, optimizer = optimizer)

    return λ, P, RG_ball_star, Δ₁x, Iₙ
end
