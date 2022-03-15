import PropertyT_new.MOI

function sos_from_matrix(Q::AbstractMatrix, support, RG::StarAlgebra)
    mn = LinearAlgebra.checksquare(Q)

    # Changing Q to the corresponding interval-entry matrix
    Q_interval = map(x->@interval(x), Symmetric((Q.+Q')./2))
    P_interval_RG = RG.(Q_interval'*Q_interval)

    m = length(support)
    n,r = divrem(mn, m)
    @assert iszero(r)
    Iₙ = [(i == j ? one(RG) : zero(RG)) for i in 1:n, j in 1:n]

    x = reshape([RG(s) for s in support], m, 1)
    xx = collect(Iₙ ⊗ x)
    result = xx' * P_interval_RG * xx

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

    residual = eoi - sos_from_matrix(Q, support, RG)
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
    λₐₚ, Pₐₚ, termination_status, RG, Δ₁, Iₙ = spectral_gaps_approximated(
        h,
        relations,
        half_basis;
        optimizer = optimizer,
    )
    @info "Termination status: " termination_status

    termination_status != MOI.OPTIMAL && return termination_status, @interval(-1)

    Qₐₚ = real(sqrt(Symmetric((Pₐₚ .+ Pₐₚ') ./ 2)))

    @info "Approximated λ: " λₐₚ

    certified_sgap = certify_sos_decomposition(Δ₁, Iₙ, λₐₚ, Qₐₚ, half_basis, RG)

    @info "Certified λ (interval atithmetic): " certified_sgap

    return termination_status, certified_sgap
end
