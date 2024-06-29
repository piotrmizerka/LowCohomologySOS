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
)
    λ_interval = @interval(λ)
    eoi = _eoi(X, λ_interval, order_unit)

    residual = eoi - sos_from_matrix(parent(first(X)), Q, support)
    l1_norm = sum(x -> norm(x, 1), residual)

    @info "l₁ norm of the error in interval arithmetic:" l1_norm radius(l1_norm)

    result = λ_interval - l1_norm

    return result.lo > 0, result
end

function spectral_gaps_certification(
    h,
    relations,
    half_basis,
    S = gens(parent(first(relations)));
    optimizer,
)
    solution = spectral_gaps_approximated(
        h,
        relations,
        half_basis,
        S;
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
    )

    @info "Certified λ (interval atithmetic): " certified_sgap

    return solution.termination_status, certified_sgap
end

# Additional function allowing to view the supporting summand factors of soses.
# Not meeded for major computations - for now just to view the pattern hopefully.
function sos_summand_factors(RG::StarAlgebra, Q::AbstractMatrix, support, cutoff, no_digits)
    result = []
    mn = LinearAlgebra.checksquare(Q)
    m = length(support)
    n = div(mn,m)
    x = kron(Matrix(I, n, n), [RG(s) for s in support])
    for i in 1:m
        Qi = [Q[(i-1)*n+k,s] for k in 1:n, s in 1:mn]
        summand_factor = Qi * x
        summand_factor_red = [0.0*zero(RG) for k in 1:n, l in 1:n]
        for k in 1:n
            for l in 1:n
                for s in 1:m
                    if abs(summand_factor[k,l](support[s])) > cutoff
                        coeff = round(summand_factor[k,l](support[s]),digits = no_digits)
                        summand_factor_red[k,l] += coeff*RG(support[s])
                    end
                end
            end
        end
        push!(result,summand_factor_red)
    end
    return result
end