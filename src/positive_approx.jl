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

function sos_problem(
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

# We want to find α, β, and λ such that
# α*A+β*B-M-λ*order_unit = SOS
function sos_problem(
    A::AbstractMatrix{<:AlgebraElement},
    B::AbstractMatrix{<:AlgebraElement},
    M::AbstractMatrix{<:AlgebraElement},
    order_unit::AbstractMatrix{<:AlgebraElement},
    upper_bound::Float64 = Inf,
)
    @assert size(A) == size(B) == size(M) == size(order_unit)
    @assert !isempty(A) && !isempty(B) && !isempty(M) && !isempty(order_unit)

    RG = parent(first(A))
    @assert all(x -> parent(x) === RG, A)
    @assert all(x -> parent(x) === RG, B)
    @assert all(x -> parent(x) === RG, M)
    @assert all(x -> parent(x) === RG, order_unit)

    m = size(RG.mstructure, 1)
    n = LinearAlgebra.checksquare(A)
    mn = m * n
    result = JuMP.Model()

    JuMP.@variable(result, P[1:mn, 1:mn], Symmetric)
    JuMP.@constraint(result, sdp, P in JuMP.PSDCone())

    JuMP.@variable(result, λ)
    JuMP.@objective(result, Max, λ)

    if upper_bound < Inf
        λ = JuMP.@constraint(result, λ <= upper_bound)
    end

    JuMP.@variable(result, α)
    JuMP.@variable(result, β)
    JuMP.@constraint(result, α >= 0.3)
    JuMP.@constraint(result, β >= 0.05)

    cnstrs = constraints(RG.mstructure)
    @assert length(cnstrs) == length(basis(RG))

    for idx in CartesianIndices(A)
        i, j = Tuple(idx)
        Pⁱʲ = @view P[KroneckerDelta{n}(i, j)]

        aij = A[idx]
        bij = B[idx]
        mij = M[idx]
        uij = order_unit[idx]

        for (RG_g, g) in zip(cnstrs, basis(RG))
            JuMP.@constraint(result, α*aij(g) + β*bij(g) - mij(g) - λ * uij(g) == dot(RG_g, Pⁱʲ))
        end
    end
    return result
end

# h:Free group --> our group G
function spectral_gap_elements(
    h,
    relations,
    half_basis,
    S = gens(parent(first(relations)))
)
    @assert !isempty(relations)
    
    G = parent(h(first(relations))) # target of h

    d₁ = jacobian_matrix(relations, S)

    Δ₁, Δ₁⁺, Δ₁⁻ = let RG = group_ring(G, half_basis, star_multiplication = false)
        d₁x = embed.(Ref(h), d₁, Ref(RG))
        d₀x = embed.(Ref(h), d₀(parent(first(d₁)), S), Ref(RG))

        Δ₁⁺ = d₁x' * d₁x
        Δ₁⁻ = d₀x * d₀x'
        Δ₁⁺ + Δ₁⁻, Δ₁⁺, Δ₁⁻
    end

    RG = group_ring(G, half_basis, star_multiplication = true)

    Δ₁x = embed.(identity, Δ₁, Ref(RG))
    Δ₁⁺x = embed.(identity, Δ₁⁺, Ref(RG))
    Δ₁⁻x = embed.(identity, Δ₁⁻, Ref(RG))

    n = length(S)
    @assert size(Δ₁x, 1) === size(Δ₁x, 2) === n
    Iₙ = [i ≠ j ? zero(RG) : one(RG) for i in 1:n, j in 1:n]

    return Δ₁x, Iₙ, Δ₁⁺x, Δ₁⁻x
end

function get_solution(m::JuMP.Model)
    λ = JuMP.value(m[:λ])
    Q = let P = JuMP.value.(m[:P])
        if any(isnan, P) || any(isinf, P)
            @error "obtained solution contains NaNs or ±Inf"
            P
        else
            real(sqrt(Symmetric((P + P')./2)))
        end
    end
    return λ, Q
end

function spectral_gaps_approximated(
    h,
    relations::AbstractVector{<:FPGroupElement},
    half_basis,
    S = gens(parent(first(relations)));
    optimizer,
)
    Δ₁x, Iₙ, Δ₁⁺x, Δ₁⁻x = spectral_gap_elements(h, relations, half_basis, S)

    Δ₁_sos_problem = sos_problem(Δ₁x, Iₙ)

    JuMP.set_optimizer(Δ₁_sos_problem, optimizer)
    JuMP.optimize!(Δ₁_sos_problem)

    status = JuMP.termination_status(Δ₁_sos_problem)

    λ, Q = get_solution(Δ₁_sos_problem)

    solution = (
        laplacian = Δ₁x,
        unit = Iₙ,
        termination_status = status,
        λ = λ,
        Q = Q
    )

    return solution
end
