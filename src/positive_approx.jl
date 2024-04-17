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
        JuMP.@constraint(result, λ <= upper_bound)
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

# h:Free group --> our group G
function spectral_gap_elements(
    h,
    relations,
    half_basis;
    S = gens(parent(first(relations))),
    twist_coeffs = true
)
    @assert !isempty(relations)
    
    G = parent(h(first(relations))) # target of h

    d₁ = jacobian_matrix(relations, S)

    Δ₁, Δ₁⁺, Δ₁⁻, d₀_ = let RG = group_ring(G, half_basis, star_multiplication = false)
        d₁x = embed.(Ref(h), d₁, Ref(RG))
        d₀x = embed.(Ref(h), d₀(parent(first(d₁)), S), Ref(RG))

        Δ₁⁺ = d₁x' * d₁x
        Δ₁⁻ = d₀x * d₀x'
        Δ₁⁺ + Δ₁⁻, Δ₁⁺, Δ₁⁻, d₀x
    end

    RG = group_ring(G, half_basis, star_multiplication = true)

    Δ₁x = twist_coeffs ? embed.(identity, Δ₁, Ref(RG)) : Δ₁
    Δ₁⁺x = twist_coeffs ? embed.(identity, Δ₁⁺, Ref(RG)) : Δ₁⁺
    Δ₁⁻x = twist_coeffs ? embed.(identity, Δ₁⁻, Ref(RG)) : Δ₁⁻

    n = length(S)
    @assert size(Δ₁x, 1) === size(Δ₁x, 2) === n
    Iₙ = [i ≠ j ? zero(parent(first(Δ₁x))) : one(parent(first(Δ₁x))) for i in 1:n, j in 1:n]

    if twist_coeffs
        return Δ₁x, Iₙ, Δ₁⁺x, Δ₁⁻x
    else
        return Δ₁x, Iₙ, Δ₁⁺x, Δ₁⁻x, d₀_
    end
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
