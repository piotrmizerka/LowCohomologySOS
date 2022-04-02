function sos_from_matrix(RG::StarAlgebra, Q::AbstractMatrix, support)
    mn = LinearAlgebra.checksquare(Q)

    # Changing Q to the corresponding interval-entry matrix
    Q_interval = map(x->@interval(x), Symmetric((Q.+Q')./2))
    P_interval = Q_interval'*Q_interval

    m = length(support)
    n,r = divrem(mn, m)
    @assert iszero(r)

    x = kron(Matrix(I, n, n), [RG(s) for s in support])
    result = permutedims(x) * P_interval * x

    return result
end

_eoi(X, λ::Number, u) = X-λ*u
_eoi(X::AbstractMatrix, λ::Number, u::AbstractMatrix) = X-Ref(λ).*u
function certify_sos_decomposition(
    X,
    order_unit,
    λ::Number,
    Q::AbstractMatrix,
    support,
    RG::StarAlgebra,
)
    λ_interval = @interval(λ)
    eoi = _eoi(X, λ_interval, order_unit)

    residual = eoi - sos_from_matrix(RG, Q, support)
    l1_norm = sum(x -> norm(x, 1), residual)

    @info "l₁ norm of the error in interval arithmetic:" l1_norm radius(l1_norm)

    result = λ_interval - l1_norm

    return result
end

function spectral_gaps_certification(
    h,
    relations,
    half_basis;
    optimizer,
)
    solution = spectral_gaps_approximated(
        h,
        relations,
        half_basis;
        optimizer = optimizer,
    )
    @info "Termination status: " solution.termination_status

    @info "Approximated λ: " solution.λ

    certified_sgap = certify_sos_decomposition(
        solution.laplacian,
        solution.unit,
        solution.λ,
        solution.Q,
        half_basis,
        parent(first(solution.laplacian)),
    )

    @info "Certified λ (interval atithmetic): " certified_sgap

    return solution.termination_status, certified_sgap
end
